#!/usr/bin/env python3
"""Reproduce the README's NVE drift and isolate its cutoff discontinuity.

Run from any directory; no third-party Python packages are required.
Only writes diagnostic inputs/results to --output. Does not change LAMMPS.
"""

import argparse
import csv
import json
import os
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]


def run_case(binary, output, name, script):
    inp = output / (name + ".in")
    inp.write_text(script)
    result = subprocess.run(
        [str(binary), "-log", "none", "-in", str(inp)],
        cwd=output, text=True, capture_output=True,
        env={**os.environ, "OMP_NUM_THREADS": "1"}, check=False,
    )
    (output / (name + ".log")).write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f"{name} failed: {result.stdout[-3000:]} {result.stderr}")
    rows = []
    active = False
    for line in result.stdout.splitlines():
        fields = line.split()
        if fields and fields[0] == "Step":
            active = True
            continue
        if line.startswith("Loop time"):
            active = False
        if active and len(fields) == 4:
            try:
                rows.append([float(x) for x in fields])
            except ValueError:
                pass
    if len(rows) < 2:
        raise RuntimeError(f"No energy trajectory in {name}")
    with (output / (name + ".csv")).open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(["step", "ke", "pe", "etotal"])
        writer.writerows(rows)
    energies = [r[3] for r in rows]
    summary = {
        "case": name,
        "initial_energy": energies[0],
        "final_energy": energies[-1],
        "delta_energy": energies[-1] - energies[0],
        "max_abs_error": max(abs(e - energies[0]) for e in energies),
        "energy_span": max(energies) - min(energies),
        "samples": len(rows),
    }
    print(json.dumps(summary), flush=True)
    return summary


def many_body(shift):
    # Same grid, seed, coefficients, alpha, and fixed z* as the earlier ec_new.in.
    # The original placed thermo_modify before thermo_style, which reset norm to
    # the LJ default (per atom). Specify that normalization explicitly here.
    return f"""units lj
atom_style atomic
boundary p p f
atom_modify map array
processors * * 1
pair_style lj/cut/mlo 9.34 0.0 2.5 5.0 theta_check off zstar fixed
read_data "{ROOT / 'mlo_tests/grid.data'}"
pair_coeff 1 1 1.7 1.00
pair_coeff 2 2 0.9 0.80
pair_coeff 1 2 1.2 0.90
pair_modify shift {shift}
group ga type 1
group gb type 2
velocity all create 0.5 9871 loop geom
timestep 1.0e-4
fix 1 all nve
fix 3 ga free_energy 1.0 3.06 0.8
fix 4 gb free_energy 1.0 10.0 0.0
fix_modify 3 energy yes
fix_modify 4 energy yes
thermo_style custom step ke pe etotal
thermo_modify norm yes format float %.17g
thermo 1
run 20000
"""


def crossing(shift, dt):
    # Saturated Theta makes the reaction-coordinate forces negligible. The pair
    # moves outward through r_c = 2.5, with no thermostat or landscape fix.
    return f"""units lj
atom_style atomic
boundary p p f
atom_modify map array
processors * * 1
region box block -10 10 -10 10 -10 10
create_box 1 box
mass 1 1.0
create_atoms 1 single 0 0 1
create_atoms 1 single 2.49993 0 1
pair_style lj/cut/mlo 4.0 0.0 2.5 50.0 zstar fixed theta_check off
pair_coeff * * 1.0 1.0
pair_modify shift {shift}
group left id 1
group right id 2
velocity left set -0.5 0 0
velocity right set 0.5 0 0
fix integrate all nve
neighbor 0.3 bin
neigh_modify every 1 delay 0 check no
timestep {dt:.17g}
thermo_style custom step ke pe etotal
thermo_modify norm no format float %.17g
thermo 1
run {round(0.001 / dt)}
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lmp", type=Path, default=ROOT / "build-new/lmp")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    binary = args.lmp.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    summaries = []
    for shift in ("no", "yes"):
        summaries.append(run_case(binary, output, f"many_body_shift_{shift}", many_body(shift)))
    for dt in (1e-4, 1e-5, 1e-6):
        for shift in ("no", "yes"):
            name = f"crossing_shift_{shift}_dt_{dt:.0e}"
            summaries.append(run_case(binary, output, name, crossing(shift, dt)))
    (output / "summary.json").write_text(json.dumps(summaries, indent=2) + "\n")
    # This checks a predicted physical failure mode, independently of the force
    # implementation: unshifted energy loses continuity at an attractive cutoff.
    expected_jump = -4.0 * (2.5 ** -12 - 2.5 ** -6)
    for item in summaries[2:]:
        if "shift_no" in item["case"]:
            assert abs(item["delta_energy"] - expected_jump) < 5e-6, item
        else:
            assert item["max_abs_error"] < 5e-6, item
    assert summaries[0]["max_abs_error"] > 1e-3, summaries[0]
    assert summaries[1]["max_abs_error"] < 1e-5, summaries[1]
    print("PASS: reproduced cutoff jump and verified its removal by shift yes.")


if __name__ == "__main__":
    main()

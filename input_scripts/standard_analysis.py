# Analyze a production trajectory in the MDA conda environment.
# LAMMPS can export one argument per line to a .analysis.args file:
# python standard_analysis.py -t trajectory.lammpstrj @run.analysis.args
# Arguments can also be supplied directly on the command line.

import argparse
import math
import pathlib

argument_parser = argparse.ArgumentParser(
    description="Standard analysis of production trajectory", fromfile_prefix_chars="@"
)
argument_parser.add_argument("--trajectory", "-t", required=True, help="Path to the trajectory file")
argument_parser.add_argument("--gamma", "-g", type=float, required=True, help="langevinMLO gamma_z/gamma_xy ratio")
argument_parser.add_argument("--epsilon", "-e", type=float, required=True, help="epsilon parameter used in this simulation")
argument_parser.add_argument("--timestep", type=float, required=True, help="integration timestep used in this simulation")
argument_parser.add_argument("--saving_interval", type=int, default=1000, help="Interval at which frames were saved in the trajectory")

argument_parser.add_argument("--barrier_left", type=float, required=True, help="positive distance from zero to the inactive basin boundary")
argument_parser.add_argument("--barrier_right", type=float, required=True, help="positive distance from zero to the active basin boundary")

arguments = argument_parser.parse_args()
for name in ("gamma", "timestep", "barrier_left", "barrier_right"):
    value = getattr(arguments, name)
    if not math.isfinite(value) or value <= 0:
        argument_parser.error(f"--{name} must be finite and positive")
if not math.isfinite(arguments.epsilon) or arguments.epsilon < 0:
    argument_parser.error("--epsilon must be finite and nonnegative")
if arguments.saving_interval <= 0:
    argument_parser.error("--saving_interval must be positive")

trajectory_path = pathlib.Path(arguments.trajectory)

if not trajectory_path.is_file():
    raise FileNotFoundError(f"Trajectory file not found: {trajectory_path}")

gamma = arguments.gamma
epsilon = arguments.epsilon
timestep = arguments.timestep
saving_interval = arguments.saving_interval
barrier_left = arguments.barrier_left
barrier_right = arguments.barrier_right

# 1) Extract the PMF for the z-coordinate

import MDAnalysis as mda
import numpy as np
import matplotlib

matplotlib.use("Agg")  # Slurm jobs have no display.
import matplotlib.pyplot as plt

u = mda.Universe(str(trajectory_path), format="LAMMPSDUMP", dt=timestep)
n_atoms = len(u.atoms)
n_frames = len(u.trajectory)
if n_atoms == 0 or n_frames < 2:
    argument_parser.error("The trajectory must contain atoms and at least two saved frames")
box_dimensions = u.dimensions.copy()


print('-' * 50)
print(f'Number of atoms: {n_atoms}')
print(f'Number of frames: {n_frames}')
print(f'Box dimensions: {box_dimensions}')
print('-' * 50)

z_positions = []
frame_steps = []
for ts in u.trajectory:
    frame_steps.append(ts.data["step"])
    for atom in u.atoms:
        z_positions.append(atom.position[2])

# Check the exported dump interval against the actual saved LAMMPS steps.
if any(second - first != saving_interval for first, second in zip(frame_steps, frame_steps[1:])):
    argument_parser.error("--saving_interval does not match the trajectory step spacing")
t = timestep * (frame_steps[-1] - frame_steps[0])

# MDAnalysis shifts the box origin to zero; these inputs use symmetric z bounds.
zdim = box_dimensions[2]
zdim_shift = zdim / 2.0
print(f"z-dimension of the box: {zdim}, shifting z-positions by {zdim_shift} to center around 0")

z_positions = np.array(z_positions) - zdim_shift

# range from -2.5 to 2.5 is plenty enough for standard phi(z) potentials

RANGE_UPPER = 2.5
RANGE_LOWER = -2.5
counts, bins = np.histogram(z_positions, bins=100, density=True, range=(RANGE_LOWER, RANGE_UPPER))
with np.errstate(divide="ignore"):
    free_energy = -np.log(counts)  # Empty bins have infinite free energy.
fig, ax = plt.subplots(1,2, figsize=(12, 5))
bin_centers = (bins[:-1] + bins[1:]) / 2
ax[0].stairs(counts, bins)
ax[1].plot(bin_centers, free_energy,lw=0, ms=5, marker='o')

ax[0].set_ylabel(r"$\mathcal{P}(z)$", fontsize=15)
ax[0].set_xlabel(r"$z$", fontsize=15)
ax[1].set_xlabel(r"$z$", fontsize=15)
ax[1].set_ylabel(r"$-kT \log[\mathcal{P}(z)]$", fontsize=15)

plt.tight_layout()

out_dir = trajectory_path.parent / "analysis" / trajectory_path.stem
out_dir.mkdir(parents=True, exist_ok=True)

# save the plots
plt.savefig(out_dir / "z_distribution.png")
plt.close(fig)

# save the histogram and the free energy profile
np.save(out_dir / "z_histogram.npy", counts)
np.save(out_dir / "z_bins.npy", bins)
np.save(out_dir / "free_energy_profile.npy", free_energy)
np.save(out_dir / "free_energy_bins.npy", bin_centers)

# 2) jumpy analysis

def assign_basin(all_atoms, z_box_dimension: float, barrier_left: float, barrier_right: float) -> tuple[mda.core.groups.AtomGroup, int, int, int, list[tuple[int, str]]]:


    barrier_z_active = (z_box_dimension / 2) + barrier_right
    barrier_z_inactive = (z_box_dimension / 2) - barrier_left

    counts = 0
    counts_active_to_inactive = 0
    counts_inactive_to_active = 0
    jump_events = []
    for atom in all_atoms:
        z = atom.position[2]
        if atom.element == 'U':
            # print(z)
            if z < barrier_z_inactive:
                atom.element = 'INACTIVE'
            elif z > barrier_z_active:
                atom.element = 'ACTIVE'
            else:
                pass
        
        elif atom.element == 'INACTIVE':
            if z > barrier_z_active:
                atom.element = 'ACTIVE'
                counts += 1
                print(f"Transition INACTIVE->ACTIVE at z={z:.2f} @ atom id {atom.id}")
                counts_inactive_to_active += 1
                jump_events.append((atom.id, 'INACTIVE->ACTIVE'))
            else:
                pass
        elif atom.element == 'ACTIVE':
            if z < barrier_z_inactive:
                atom.element = 'INACTIVE'
                counts += 1
                print(f"Transition ACTIVE->INACTIVE at z={z:.2f} @ atom id {atom.id}")
                counts_active_to_inactive += 1
                jump_events.append((atom.id, 'ACTIVE->INACTIVE'))
            else:
                pass
        else:
            raise ValueError(f"Unknown basin element: {atom.element}")

    assert counts == counts_active_to_inactive + counts_inactive_to_active, "Total counts should equal sum of transition counts"

    return (all_atoms, counts, counts_active_to_inactive, counts_inactive_to_active, jump_events)


# Initial assignment

initial_values = np.array(['U'] * n_atoms, dtype=object)
u.add_TopologyAttr('element', initial_values)
selection = u.select_atoms('all')

f_a = []
f_b = []
f_u = []
counts = 0
counts_active_to_inactive = 0
counts_inactive_to_active = 0
# The first frame assigns initial states without counting a transition.
for ts in u.trajectory:
    selection, c, c_ati, c_ita, jump_events = assign_basin(
        selection, zdim, barrier_left, barrier_right
    )
    counts += c
    counts_active_to_inactive += c_ati
    counts_inactive_to_active += c_ita

    group_active = u.select_atoms('element ACTIVE')
    group_inactive = u.select_atoms('element INACTIVE')
    group_unassigned = u.select_atoms('element U')

    f_a.append(group_active.n_atoms / n_atoms)
    f_b.append(group_inactive.n_atoms / n_atoms)
    f_u.append(group_unassigned.n_atoms / n_atoms)


print(counts)
print(counts_active_to_inactive)
print(counts_inactive_to_active)

print(f"dt = {timestep}")
print(f"Jumping Rate J = {counts/t / len(u.atoms)}")
print(f"Average fraction active: {np.average(f_a)}")
print(f"Average fraction inactive: {np.average(f_b)}")
print(f"Average fraction unassigned: {np.average(f_u)}")

assert counts == counts_active_to_inactive + counts_inactive_to_active, "Total counts should equal sum of transition counts"

# Save the summary printed above in a human-readable format as well.
results_path = out_dir / "analysis_results.txt"
with results_path.open("w", encoding="utf-8") as results_file:
    results_file.write("Standard analysis results\n")
    results_file.write("=" * 26 + "\n")
    results_file.write(f"Trajectory: {trajectory_path}\n")
    results_file.write(f"Number of atoms: {n_atoms}\n")
    results_file.write(f"Number of frames: {n_frames}\n")
    results_file.write(f"Gamma: {gamma}\n")
    results_file.write(f"Epsilon: {epsilon}\n")
    results_file.write(f"Saving interval: {saving_interval}\n")
    results_file.write(f"Inactive basin: z < {-barrier_left}\n")
    results_file.write(f"Active basin: z > {barrier_right}\n")
    results_file.write(f"Box dimensions: {box_dimensions}\n")
    results_file.write(f"z-dimension of the box: {zdim}\n")
    results_file.write(f"Total transitions: {counts}\n")
    results_file.write(f"ACTIVE -> INACTIVE transitions: {counts_active_to_inactive}\n")
    results_file.write(f"INACTIVE -> ACTIVE transitions: {counts_inactive_to_active}\n")
    results_file.write(f"Timestep: {timestep}\n")
    results_file.write(f"Total simulation time: {t}\n")
    results_file.write(f"Jumping Rate J: {counts / t / len(u.atoms)}\n")
    results_file.write(f"Average fraction active: {np.average(f_a)}\n")
    results_file.write(f"Average fraction inactive: {np.average(f_b)}\n")
    results_file.write(f"Average fraction unassigned: {np.average(f_u)}\n")

print(f"Analysis saved to: {out_dir}")


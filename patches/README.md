## MLO Changes 

### Added
- (2026/04) added `fix free_energy` for the constant potential
- (2026/01) Added anisotropic Langevin thermostat `langevinMLO` where $\gamma_{xy} \neq \gamma_z$
- (2026/01) Added modified AH potential `lj/cut/mlo` where $\lambda = \lambda(z)$ and the LJ part only depends on the distance in the $xy$-plane.

## Bugfixes

- (2026/04) $\lambda_{ij}$ is now correctly calculated as $\lambda_{ij} = \Theta_i \Theta_j$
- (2026/04) self-interactions

### 2026/08 — code review pass

Fixes for the findings in [`CODE_REVIEW.md`](../../CODE_REVIEW.md) (dated 2026-08-31).
Every change is marked in the source with a `FIXED by CLAUDE (2026-08-31)` comment.
`patches/` and `src/` are kept byte-identical.

**⚠️ Two of these change results. See [Action items](#action-items-2026-08).**

#### `pair_lj_cut_mlo.cpp` (791 → 595 lines)

| Sev | Fix |
|---|---|
| **P1** | **The in-plane cutoff is now actually applied.** The gate tested the 3D `cutsq` (= `cut_global²`) against a distance whose z component had already been zeroed, so the effective xy interaction range was `cut_global` = 5.59 $\sigma$ instead of the declared `cutoff_2d` = 2.5 $\sigma$; `cutoff_2d` was a parsed-but-dead parameter. `cut_global` remains the *neighbor* cutoff, which stays inflated on purpose. |
| **P1** | **Uninitialized reads with more than one atom type** (not in the review). `init_one` mirrored `lj1`–`lj4` and `offset` into the lower triangle but not `epsilon`, `sigma`, or `LJ_MINIMUM` — which `compute()` indexes as `[itype][jtype]` with no ordering guarantee. Every pair with `itype > jtype` read uninitialized memory for the well depth and the AH crossover radius. Affects `two_particles.in` (2 types). |
| P2 | **Ashbaugh–Hatch shift made consistent.** Both branches are now $U = U_{LJ} + (1-\lambda)\varepsilon - \lambda\,U_{LJ}(r_c)$ for $r \le 2^{1/6}\sigma$ and $U = \lambda\,(U_{LJ} - U_{LJ}(r_c))$ above it, with the z-force as $-(\partial U/\partial\lambda)\,\partial\lambda/\partial z$ in both. Previously the repulsive branch carried no shift, and the attractive branch shifted its energy by `offset` but its force by `λ·offset`, so the force was not the gradient of the energy. `offset` is now computed from `cutoff_2d`, not `cut[i][j]`. |
| P2 | $\alpha$ (steepness of $\Theta(z)$) is now an **optional 4th `pair_style` argument**, default 50 — it used to be hardcoded, and had drifted to 100 in the `single()` code path. |
| P2 | **Geometric assumptions are checked** in a new `init_style()`: (a) $\sqrt{c_{2d}^2 + \Delta z_{max}^2} \le c_{global}$, so no in-plane neighbor is culled by its z separation; (b) with periodic z, $L_z > 2(\text{cutoff}+\text{skin})$, so no particle meets its own z image at zero in-plane distance (infinite force). Both were previously load-bearing but unwritten. |
| P2 | `single()` returned a hardcoded `10.0` with `single_enable = 1`. Now `single_enable = 0` and an explicit error — `fix widom`/`fix gcmc` auto-promote to the full-energy path. `born_matrix()` had the same defect (it evaluated the plain 3D LJ derivatives) and is disabled the same way. |
| P3 | `1.0/rsq_2d` was evaluated for every neighbor before any cutoff test. Moved inside the gate (a real speedup now that most neighbors are culled) with a zero-distance guard. |
| P3 | Removed ~200 lines of dead debug scaffolding: the `pure_lj_interactions` and `output_mlo_interactions` blocks, the debug counters, the two unreachable assertion `printf`s, the version banner in `init_one`, and `cut_respa`. |
| P3 | `pair_modify tail yes` is now rejected — the 3D tail-correction formula does not apply to an xy-only, $\lambda$-scaled interaction. |
| — | `thermo press` remains meaningless (the global virial comes from `virial_fdotr_compute()`, which also picks up the non-pairwise $\lambda(z)$ z-forces, and the box volume contains $L_z$). Not fixable inside the pair style; documented at the `ev_tally` call. |

*Verification:* the rewritten potential was checked against finite differences — max relative
gradient error **2.7e-10** with the shift on and off, energy continuous at the AH crossover,
and $U(r_c) \approx 5\times10^{-18}$ for arbitrary $\lambda$ under `pair_modify shift yes`
(the old form did not vanish at the cutoff).

#### `fix_free_energy.cpp`

| Sev | Fix |
|---|---|
| P2 | `post_force` **no longer mutates velocities**. `v_z` is pinned once in `setup()`/`min_setup()` and re-pinned in a new `post_integrate()`. `FixWidom::energy_full()` calls `modify->post_force()` on every insertion trial, which used to corrupt the real trajectory's z-velocities under `disable_reactions`. |
| P2 | `disable_reactions` **no longer depends silently on fix definition order**. `init()` walks the fix list and errors out, naming the offending fix, if any later fix carries `POST_FORCE`. |
| P2 | `MIN_POST_FORCE` added, so $\phi(z)$ is **active during `minimize`**. Without it the pair style's z-force had nothing to oppose it and minimization drove every particle to $\lambda = 1$. Confirmed: with the fix defined before `minimize`, $\phi$ per atom settles at exactly $-b^2/4 = -2.3409$ (every atom in a well). |
| P3 | `compute_scalar()` reset its own cache flag on every call, so each call did a fresh `MPI_Allreduce`. `eflag` is now cleared in `post_force` instead. |
| P3 | Implemented `dof()`, so temperature computes know about the frozen z degree of freedom under `disable_reactions` (`compute temp` was reading ⅔ of the true in-plane temperature). **Note:** this makes plain `compute temp` correct, but `compute temp/partial 1 1 0` then over-reports by 3/2, because `ComputeTempPartial` prorates fix-removed DOF evenly across active dimensions — right for SHAKE-like constraints, wrong for a purely-z freeze. Use plain `compute temp` when `disable_reactions` is on. |
| P3 | `pow(z,3)`/`pow(z,4)` replaced by multiplies; comment added that `a` multiplies only the quartic term and is **not** an overall scale of $\phi(z)$, despite the input scripts naming it `s`. |

#### `fix_langevinmlo.cpp`

| Sev | Fix |
|---|---|
| P2 | **Argument check was off by one:** `narg < 7` guarded a syntax that reads `arg[7]`. Now `narg < 8`, plus a keyword loop that *rejects* unrecognized extra arguments instead of ignoring them. The seed's error-annotation index was corrected from 6 to 7. |
| P3 | **`tally yes/no` parsing restored.** `tallyflag` was permanently 0, so `END_OF_STEP` was never masked and `compute_scalar()` always returned 0 while `ecouple_flag = 1` advertised otherwise. `ecouple_flag` is now set to `tallyflag`, so `Ecouple`/`Econserve` either report the real number or don't count this fix at all. |
| P3 | Removed the dead `franprev`/`lv` machinery — `grow_arrays`, `copy_arrays`, `pack_exchange`, `unpack_exchange` and both arrays, which were never allocated, never registered via an atom callback, and never freed. `memory_usage()` now reports only arrays that exist. |
| P3 | Dropped two unused leftover locals in `post_force_templated` (the only warnings in the file, multiplied across all 32 template instantiations). |

`tally yes` together with `thermo_style custom … econserve` is the direct diagnostic for the
timestep problem below.

## Changes

### Interface changes (2026/08)

```
pair_style lj/cut/mlo <neighbor cutoff> <z_star> <xy cutoff> [alpha]
fix ID group langevinMLO Tstart Tstop damp aniso seed [tally yes/no]
```

Both new arguments are optional and default to the previous behaviour
(`alpha = 50`, `tally no`), so all five existing input scripts run unchanged.

`z_star`, `cutoff_2d` and `alpha` are **not** stored in restart files. Adding them would
break the existing `after_equil.out` / `state_after_probing_chem_reactions.out`, so the
format is unchanged and `init_style()` errors out loudly if no `pair_style` command was
issued after `read_restart`. The input scripts already reissue it.

> **Superseded by 2026/09.** These settings *are* stored now, behind a versioned header that
> still reads the old files. See [Restart format](#restart-format-2026-09).

<a name="action-items-2026-08"></a>
### Action items

1. **Re-run the phase behaviour (Figure 3).** Every run to date used an attractive range of
   5.59 $\sigma$ in the plane rather than the 2.5 $\sigma$ the input scripts declare. That
   changes the second virial coefficient, $T_c$, and the coexistence curve. To reproduce the
   old behaviour for comparison, pass the xy cutoff equal to the neighbor cutoff.
2. **Re-check anything from `two_particles.in`** — it has two atom types and was reading
   uninitialized memory for `epsilon`/`LJ_MINIMUM` on roughly half its pairs.
3. **`test.in` and `widom.in` use a timestep ~10× too large** for $\Gamma = 100$:
   $\Delta t\,\gamma_z/m = 0.77$, where the explicit Langevin kick needs $\lesssim 0.1$.
   `two_particles.in` and `probe_chem_reactions.in` are fine at 0.077. Not a code bug — fix
   in the input scripts. Validate with `compute temp/partial 0 0 1` vs `1 1 0` and by
   checking $p(z)$ against $e^{-\beta\phi(z)}$.
4. **Fix the $\Gamma$ convention** in the progress report and `PROJECT_BASELINE.md`. The code
   implements $\gamma_z = \Gamma\gamma_{xy}$, hence $\Gamma = D_{xy}/D_z$ — the inverse of
   the report's stated definition. `GAMMA = 100` means $D_z$ is 100× *smaller*.
5. Still open, from the review: whether $\lambda_{ij} = \Theta_i\Theta_j$ (product) or
   $(\Theta_i+\Theta_j)/2$ (standard AH) is intended; whether `z_star` should be solved from
   $4z^3 - 2bz + f = 0$ rather than supplied by hand; and whether the zero $\varepsilon$ for
   the fuel species in `two_particles.in` is a placeholder.

**Status after 2026/09:** items 2 and 4 are still open. Item 1 is unchanged (the in-plane
cutoff is now per type pair, but the *value* the old runs used was still 5.59 $\sigma$).
Item 3 is unchanged — still an input-script matter. Of item 5, two of the three are now
answered: `z_star` **is** solved from $4z^3 - 2bz + f = 0$, and the zero $\varepsilon$ was
**not** a placeholder but a bug — it deleted the fuel's excluded volume. Whether
$\lambda_{ij}$ should be the product or the arithmetic mean is still a modelling choice for
the paper.


### 2026/09 — per-type-pair interaction matrix

Goal: make `lj/cut/mlo` safe to drive with a full interaction matrix (different LJ parameters
*and* different attraction for every atom-type pair), and turn "inactive states never attract"
from an assumption into something the code checks.

`patches/` and `src/` are kept byte-identical (only `src/` is compiled — CMake globs
`src/*.cpp`, so a stale `patches/` copy is silently wrong).

**⚠️ `two_particles.in` results change by a factor of ~2, for a good reason. See
[Action items](#action-items-2026-09).**

#### The potential now

With $\Lambda_{ij}$ = `lambda_max[i][j]`, $\Theta_i = 1/(1+e^{-\alpha(z_i - z^*_{t_i})})$,
$r^2 = \Delta x^2 + \Delta y^2$ and $r_c \equiv c_{2d,ij}$ the per-type-pair in-plane cutoff:

$$\lambda_{ij} = \Lambda_{ij}\,\Theta(z_i)\,\Theta(z_j)$$

The following continuous-energy form requires **`pair_modify shift yes`** and
$r_c \ge 2^{1/6}\sigma_{ij}$, as in the current inputs. The shift is not enabled by
default; see [NVE drift resolved](#nve-drift-resolved-2026-09).

$$U = \begin{cases} U_{LJ}(r) + (1-\lambda_{ij})\varepsilon_{ij} - \lambda_{ij} U_{LJ}(r_c) & r < r_c,\ r \le 2^{1/6}\sigma_{ij} \\ \lambda_{ij}\left(U_{LJ}(r) - U_{LJ}(r_c)\right) & 2^{1/6}\sigma_{ij} < r < r_c \\ 0 & r \ge r_c\end{cases}$$

$\varepsilon_{ij}$ owns the repulsive core, $\Lambda_{ij}$ owns the attraction. So
$\Lambda_{ij} = 0$ is **exact WCA at full $\varepsilon_{ij}$** — excluded volume with no
attraction ever, which was previously inexpressible — and $\Lambda_{ij} = 1$ reproduces the
2026/08 potential bit for bit.

#### What was wrong

| Sev | Fix |
|---|---|
| **P1** | **$\varepsilon_{ij}$ controlled the core *and* the well depth**, so a repulsion-only pair could not be written down. `two_particles.in:100-101` already attempted it (`pair_coeff 1 2 0.0 0.75`, `pair_coeff 2 2 0.0 0.5`): $\varepsilon = 0$ zeroes `lj1`–`lj4` *and* `offset`, so all 64 fuel particles had **no interaction whatsoever** and could interpenetrate freely. Fixed by the new `lambda_max` matrix. |
| **P1** | **Periodic-image ghosts carry the wrong chemical state** (not in the 2026/08 review). `AtomVec::pack_comm` stores a z-image ghost at `x[j][2] + pbc[2]*zprd`, and `compute()` reads `x[j][2]` to evaluate $\Theta(z_j)$. A ghost of an atom in the inactive well is therefore seen at $z + L_z$, where $\Theta \approx 1$ instead of $\approx 0$: its chemical state, and the sign of its attraction, silently inverts. Dormant only because $\lvert z\rvert \lesssim 2.5 \ll 20 - 5.89$ — a much tighter condition than the 2026/08 $L_z > 2(\text{cutoff}+\text{skin})$ check. A wrapped $z$ is just as bad ($\Theta$ jumps, $\phi(z)$ delivers a delta-function force). Now **`boundary p p f` is required**, which removes the class; `zperiodic_ok` overrides at your own risk. |
| **P2** | **The in-plane cutoff is per type pair** (`cut_2d[i][j]`). It used to be one global scalar while `pair_coeff`'s optional 5th argument set `cut[i][j]` → `init_one` return → `cutsq[i][j]` → `cutneighsq[i][j]`, i.e. a per-pair *neighbour* cutoff that `init_style()` did not check. No script passed 5 arguments, so this was latent. A bare 5th positional is now a **hard error**: the 2026/08 `write_data_all` emitted the (inflated) neighbour cutoff in that slot, so silently reinterpreting it would reproduce the 5.59 $\sigma$ in-plane range of 2026/08 P1. |
| **P2** | **$z^*$ is per atom type and derived from the landscape.** $\Theta$'s crossover must sit at the barrier top of $\phi(z)$, and the landscape is already per type: `two_particles.in:130-131` runs two `fix free_energy` instances on type-based groups whose barrier tops are $0.1322$ and $0$, while $\Theta$ used $z^* = 0$ for both. `init_style()` now maps each atom type to its landscape (**from the atoms**, so a group defined by region or id still works), solves $\phi'(z) = 4az^3 - 2bz + f = 0$ by bisection, and puts $z^*_t$ at the barrier top. A type claimed by two landscapes is an error; a type covered by none warns (its $z$ is unbounded, because the pair z-force is $\ge 0$ identically). |
| **P2** | **The $\lvert\Delta z\rvert$ budget is now a contract, checked at every rebuild.** 2026/08 sampled the *instantaneous* z spread at setup and treated that as the requirement — backwards, since $z$ is dynamical. `z_budget` is now first class, each pair's neighbour cutoff is **derived** as $\sqrt{c_{2d,ij}^2 + z_{budget}^2}$ (so shrinking an in-plane cutoff can no longer eat the z headroom), and `compute()` rechecks the spread whenever the list is rebuilt. Cost is **zero MPI calls**: with `procgrid[2] == 1` and non-periodic z, `Comm::setup` leaves `maxneed[2] == 0`, so a rank's local+ghost set holds every possible in-plane neighbour of its locals at any $z$ and the rank-local spread is exact. |
| **P2** | **`processors * * 1` is required.** $z$ must not be domain-decomposed — it is what makes the guard above exact, and every atom sits near $z = 0$ anyway (also 2026/08 review open question 6, for load balance). |
| **P2** | **The inactive-state guarantee is now measured.** $\lambda_{ij} = \Lambda\Theta_i\Theta_j$ is never exactly zero, so the claim is about magnitude — and $\alpha$ became user-facing in 2026/08 with nothing checking it. `init_style()` reports $z_0^-, z_0^0, z_0^+$, the barrier $f^\ddagger$, $1/\alpha$, $\Theta(z_0^-)$ and the residual $\lambda_{res} = \Lambda\Theta(z_0^-)^2$ per landscape; warns above `theta_check` (default 1e-4) and errors above 100×. At $\alpha = 50$ the residual is $8\times10^{-63}$; at $\alpha = 2$ it is $2.9\times10^{-3}$ (0.3% of $\varepsilon$ per neighbour) and warns; at $\alpha = 0.5$ it errors. |
| P3 | Restart files now carry every `pair_style` setting plus `cut_2d` and `lambda_max`, behind a 64-bit magic (see [Restart format](#restart-format-2026-09)). |
| P3 | `extract()` exposes `lambda_max` (alias `lambda`, dim 2) so `fix adapt` can ramp the attraction matrix — the natural knob for the planned $\Delta\mu$ drive — plus `z_star`, `alpha`, `z_budget` (dim 0). `dim` is now set **per key**; it used to be set once to 2 at the top. Cutoffs are deliberately **not** exposed: `Pair::reinit()` calls `init_one()` but never re-runs `Neighbor::init()`, so adapting an interaction range would leave `cutneighsq` stale and silently drop pairs. |
| P3 | `init_one()` mirrors `lambda_max`, `cut_2d`, `cut_2d_sq` into the lower triangle. Required for `fix adapt`, which writes only the upper triangle before calling `reinit()` — the same class of bug as 2026/08 P1. Note `reinit()` also re-mixes any pair with `setflag == 0`, so adapting a *mixed* pair is silently a no-op: adapt the diagonals and let the mixing rule propagate. |
| P3 | `epsilon` $\ge 0$, `sigma` $> 0$, `cut_2d` $> 0$, `lambda` $\ge 0$ are validated; `lambda > 1` warns rather than errors (it is well posed — the well is simply deeper than the core amplitude — but the name implies $[0,1]$). A cutoff inside the LJ minimum warns. A zero z budget errors. |
| P3 | `fix free_energy` gained `extract()` for `a`, `b`, `f` (dim 0). No behaviour change; the pair style needs it to locate the barrier top. |
| P3 | Deleted the stale `// To-Do: Disable respa` comment — `Pair::respa_enable` defaults to 0 and `src/respa.cpp:306` already errors, so there was nothing to do. |

`compute pair/local` and `pair_write` also route through `single()`, so they error out too
(since 2026/08). Use the `full_energy` path where one exists.

#### Mixing rules

| quantity | rule when `setflag[i][j] == 0` |
|---|---|
| `epsilon` | `mix_energy` (honours `pair_modify mix`) |
| `sigma`, `cut_2d` | `mix_distance` |
| `lambda_max` | **geometric** $\sqrt{\Lambda_{ii}\Lambda_{jj}}$ by default; `lambda_mix arithmetic` for the published AH/HPS rule |

Geometric is the default deliberately. It is self-consistent with the *product* form already
chosen — $\lambda_{ij} = (\sqrt{\Lambda_i}\Theta_i)(\sqrt{\Lambda_j}\Theta_j)$, making
$\Lambda$ a per-type "bonding valence" — and it makes $\Lambda_{ii} = 0$ mean **inert with
every partner** under wildcard mixing. The arithmetic rule would give $\Lambda_{12} = 0.5$ for
an inert type paired with a fully attractive one, silently re-attracting a species declared
inert. That is exactly the failure this change exists to prevent.

#### Why $\Theta$ is *not* clamped to zero

A hard floor on $\Theta$ would make the residual identically zero, and it was rejected: it
makes $\Theta'$ discontinuous, so $U$ is no longer $C^1$, symplectic energy conservation
degrades to $O(1)$ drift per barrier crossing, and the Langevin stationary distribution is
biased. Worse, the kink would sit exactly where particles cross it $\sim k$ times per unit
time, so the error would scale with the reaction rate being measured. The analytic sigmoid
plus the diagnostic above is the answer. If an *identically*-zero residual is ever needed for
a formal claim, use a $C^2$ compactly-supported switch
$\Theta(u) = u^3(10 - 15u + 6u^2)$ on $u = (z - z^* + w)/2w$ with $w = 15/4\alpha$ (matches the
logistic slope at $z^*$) — designed but **not implemented**, since $8\times10^{-63}$ is
already far below double precision.

#### Interface changes (2026/09)

```
pair_style lj/cut/mlo <neighbor cutoff> <z_star> <xy cutoff> [alpha] [keyword values ...]

    alpha        <a>                     steepness of Theta(z)                 (default 50)
    zbudget      <dz>                    guaranteed |dz| headroom  (default sqrt(cn^2-cxy^2))
    zstar        auto | fixed            derive z* from the landscape, or use the argument
    zcheck       off | warn | strict     action of the per-rebuild |dz| guard   (default warn)
    theta_check  <tol> | off             tolerance on the residual inactive lambda (1e-4)
    lambda_mix   geometric | arithmetic  mixing for lambda_max          (default geometric)
    zperiodic_ok                         accept periodic z and its hazards

pair_coeff i j <eps> <sigma> [cut <cut_2d>] [lambda <lambda_max>]
```

The three positional arguments and the legacy 4th positional `alpha` are unchanged, so
existing `pair_style` lines still parse. `pair_coeff` defaults to `lambda 1.0` and
`cut <global xy cutoff>`, so existing 4-argument lines are unchanged too. Input scripts now
additionally need `boundary p p f` and `processors * * 1`.

<a name="restart-format-2026-09"></a>
#### Restart format

`write_restart_settings` now emits a 64-bit magic (`"MLO_PAIR"`) and a version int, then every
global; the per-pair block carries `cut_2d` and `lambda_max` instead of the derived `cut`.
`read_restart_settings` reads the magic and `fseek`s back to the old layout if it does not
match, so **`after_equil.out` and `state_after_probing_chem_reactions.out` still load** (with
a warning, and `pair_style` must still be reissued for those). A 64-bit magic rather than
32-bit because a pre-versioning file's first bytes are the double `cut_global`, whose low
mantissa word could collide.

A v2 file needs no reissued `pair_style`: `read_restart` → `run 0` reproduces `evdwl` bitwise.

#### Verification

`mlo_tests/run_tests.sh [build-dir]` — 5 checks, all passing. `mlo_tests/gen_grid.py`
regenerates the test configuration (atoms on an *xy grid*: `create_atoms random … overlap`
rejects overlaps in 3D, which does **not** guarantee in-plane separation for this pair style).

1. **$\Lambda = 0$ reproduces stock `lj/cut` WCA bitwise** at 14 radii, $U$ exactly `0.0`
   beyond $2^{1/6}\sigma$, $f_z \equiv 0$ everywhere — and the core survives ($U = +37.4$ at
   $r = 0.8\sigma$, where the old $\varepsilon = 0$ encoding gave zero force).
2. **$\Lambda = 1$, $\Theta = 1$ reproduces stock shifted `lj/cut` bitwise** at all 14 radii,
   including across the AH crossover (the branches agree to $1.5\times10^{-15}$ there).
3. **Finite-difference gradient check** via `fix numdiff` (hence `PKG_EXTRA-FIX`):
   `maxerr_z` $\le 9.5\times10^{-9}$ with $\max\lvert f_z\rvert = 1.25$. **This test must run
   at $\alpha \approx 5$.** At $\alpha = 50$ the z-forces are numerically zero away from the
   crossover and the test is vacuous — it passes even with $\Lambda$ deleted from them.
   Confirmed by injecting exactly that bug: `maxerr_z` went from $10^{-11}$ to **0.28–0.67**
   while `maxerr_y` stayed clean. Neither the regression gate ($\Lambda = 1$ is the identity)
   nor checks 1–2 ($\Theta = 1$ or $\Lambda = 0$ both zero the z-force) can catch it.
4. **Asymmetry invariant:** $\sum f_x = -8.9\times10^{-16}$, $\sum f_y \approx 0$ ($x,y$ are
   cyclic), $\sum f_z = 57.597$ **nonzero** — the model requires this, so a future "fix" that
   symmetrises the z-force fails here.
5. **MPI / newton consistency:** $\sum f_z$ = `5.7597123037901e+01` identical to 13 digits
   across `np = 1,2,4` × `newton on/off`. This is what exercises reverse communication of the
   *asymmetric* z-forces.

Also checked interactively: the landscape solver returns wells at $\pm\sqrt5$ and a barrier of
exactly $b^2/4a = 25\,kT$ for $a=1, b=10, f=0$; every guard and error path fires; the per-pair
cutoff bites (a pair at $r = 1.4$ with `cut 1.25` gives `evdwl` exactly 0 vs $-0.0979$ at
`cut 2.5`); `write_data … pair ij` round-trips bitwise; and `fix adapt` ramping `lambda_max`
reaches exactly the WCA reference.

**Regression:** a baseline binary was rebuilt from the 2026/08 `src/` (the checked-in `build/`
is stale — its binary predates the 2026/08 fixes and its CMake cache points at a deleted
directory). All four scripts are **bit-identical** across the `cut_2d` and `lambda_max` steps,
and remain bit-identical with the new code running under legacy settings
(`zstar fixed zperiodic_ok theta_check off` + `boundary p p p`). The new machinery is inert;
only the two deliberate corrections below change anything.

#### Results that change

| script | E_pair before | after | cause |
|---|---|---|---|
| `test.in` | -5.8070318 | -5.8070318 | unchanged ($b=3.06, f=0$ already had $z_0^0 = 0$) |
| `widom.in` | -1.4880694 | -1.4880694 | unchanged ($b=2, f=0$) |
| `probe_chem_reactions.in` | -1.9202685 | -1.9201233 | $z^*$ 0.256 → 0.256073 (7.6e-5) |
| `two_particles.in` | -3.3270894 | **-1.6908272** | fuel regains its excluded volume; $z^*$ 0 → 0.13223 |

#### Input script changes

- All four: `boundary p p f` and `processors * * 1`.
- `probe_chem_reactions.in`: `change_box all boundary p p f` after `read_restart`, since the
  restart file was written with a periodic z.
- `two_particles.in`: the interaction matrix now says what it meant —
  `pair_coeff 1 2 1.0 0.75 lambda 0.0` (excluded volume, never attracts) rather than
  `pair_coeff 1 2 0.0 0.75` (nothing at all). `z_max` 5.0 → 6.0: measured over 600k steps the
  z spread plateaus at **4.34**, so 5.0 left only 0.66 of headroom. And both `free_energy`
  fixes moved **above** `minimize` — the pair z-force is $\ge 0$ identically, so without
  $\phi(z)$ to oppose it the minimizer drives every particle toward $\lambda = 1$ and, at small
  $\alpha$, clean out of the z budget ($z_{max} \approx z^* + \alpha^{-1}\ln(\varepsilon\alpha/f_{tol})$).

<a name="nve-drift-resolved-2026-09"></a>
#### 2026-09-07 — NVE drift resolved

**The cause of the previously unexplained NVE drift is confirmed, and the correction
is verified.** The original NVE diagnostic and all four production inputs omitted
`pair_modify shift yes`, leaving the default `shift no`. With that setting,
`init_one()` sets `offset = 0`, while `compute()` drops the interaction at the
in-plane cutoff. Thus the attractive pair energy jumps from
$\lambda_{ij}U_{LJ}(r_c)$ just inside the cutoff to zero outside it, without a
compensating kinetic-energy change. This was an input/default mismatch with the
continuous potential written above.

**Required setting for simulations using `lj/cut/mlo`:** place the following after
the pair setup and any `read_restart`, before minimization or dynamics:

```lammps
pair_modify shift yes
```

The existing C++ implementation already applies the shift consistently in both
energy branches and the reaction-coordinate forces. No rebuild or change to
`fix nve`, `fix free_energy`, or `fix langevinMLO` is needed for this correction.

**Verification:** the preserved original case reproduced the reported
$-4.7303\times10^{-3}$ drift. With the same 196-atom initial configuration, seed,
coefficients, $\alpha = 5$, and both landscape energies included, the results over
20,000 NVE steps at $\Delta t = 10^{-4}$ were:

| Energy error per atom (LJ units) | `shift no` | `shift yes` |
|---|---:|---:|
| Final minus initial | -4.7303175944e-3 | -6.3377459836e-8 |
| Maximum absolute deviation from the initial energy | 5.7013040926e-3 | 2.3979973485e-7 |

Every step was sampled. Each error is measured against that run's own initial
energy; the maximum deviation falls by approximately **24,000 times**. A separate
two-particle cutoff-crossing test with saturated $\Theta$ confirmed the mechanism:
at $r_c = 2.5\sigma$ the unshifted outward-crossing error approaches
$0.016316891136\,\varepsilon$ as the timestep decreases from $10^{-4}$ to $10^{-6}$,
while the shifted error becomes small. This short test avoids comparing divergent
many-body trajectories at different timesteps.

The earlier finite-difference checks missed the discontinuity because they checked
forces away from the cutoff; the checked-in `numdiff.in` also explicitly enables
`shift yes`. Both landscape energies were already accounted for, and the original
NVE case had no Langevin thermostat. The suspected thermostat and strong
$\Theta'$ kick are therefore superseded by this confirmed diagnosis.

**Physical interpretation:** the shifted potential is a reasonable finite-range
interaction model: its energy approaches zero continuously at the cutoff, it
preserves the repulsive core, and attraction still depends on both chemical states.
For $r_c = 2.5\sigma$, the attractive well is about **1.63% shallower**, with depth
$0.983683\,\lambda_{ij}\varepsilon_{ij}$. Because $\lambda_{ij}$ depends on $z$, the
shift also changes the chemical forces by
$\Delta F_{z_i} = U_{LJ}(r_c)\,\partial\lambda_{ij}/\partial z_i$ inside the cutoff;
the code includes this term. Applying only a correction to recorded energies would
be insufficient. Reaction-rate results should be revalidated with the shifted
model. Small finite-timestep errors remain because the in-plane force still has an
abrupt cutoff; this is separate from the resolved large energy jump.

The [full investigation](../../NVE_DRIFT_REPORT.md),
[diagnostic script](../mlo_tests/nve_drift_diagnostic.py), and
[measured results](../mlo_tests/nve_drift_results/summary.json) are retained. To
reproduce from `lammps-stable_22Jul2025/`:

```bash
python3 mlo_tests/nve_drift_diagnostic.py --output mlo_tests/nve_drift_results
```

<a name="action-items-2026-09"></a>
#### Action items

1. **Re-run anything from `two_particles.in`.** Its 64 fuel particles previously had no
   excluded volume at all, so they were ideal-gas tracers rather than a fuel species. E_pair
   changes by a factor of ~2.
2. **Decide the $\Lambda$ mixing rule for the paper** and state the interaction law as
   $\lambda_{ij} = \Lambda^{(t_i t_j)}\Theta(z_i)\Theta(z_j)$ with that rule.
   `PROJECT_BASELINE.md` eq. (4) currently says only $\lambda \equiv \Theta(z)$.
3. **Fix the $\alpha$ convention** in the progress report and `PROJECT_BASELINE.md`. Eq. (3)
   defines $\Theta = \frac12(\tanh(\alpha(z-z_0^0))+1)$, but the code uses the logistic
   $1/(1+e^{-\alpha(z-z^*)}) = \frac12(\tanh(\alpha(z-z^*)/2)+1)$. Equal slope at $z^*$ requires
   $\alpha_{code} = 2\alpha_{report}$, so `alpha 50` in the code is $\alpha = 25$ in the
   report's convention. Fix the **report**, not the code.
4. Not done, deliberately: the `theta_form switch` indicator (above), and reverting the bare
   `delz = 0.0;` hack in stock `src/pair_lj_cut.cpp` — which is applied in `compute()` only,
   so `compute_inner/middle/outer()` and `single()` still use the full 3D `delz`, and there is
   no `rsq == 0` guard. A plain 2D LJ reference run can instead be had from
   `pair_style lj/cut/mlo <cut> -1e6 <cut_xy>`, which clamps the sigmoid to $\Theta \equiv 1$.

#### Building

The checked-in `build/` cannot be rebuilt (its cache points at a path that no longer exists).
`PKG_EXTRA-FIX` is new and load-bearing — it provides `fix numdiff`, which is the only
practical way to validate the asymmetric reaction-coordinate forces:

```
cmake -S cmake -B build-new -D CMAKE_BUILD_TYPE=RelWithDebInfo \
  -D BUILD_MPI=on -D BUILD_OMP=on -D BUILD_SHARED_LIBS=off \
  -D LAMMPS_SIZES=smallbig -D LAMMPS_MEMALIGN=64 -D PKG_MC=on -D PKG_EXTRA-FIX=on
cmake --build build-new -j
./mlo_tests/run_tests.sh build-new
```

----------------------------------------------------------------------
----------------------------------------------------------------------

This is the LAMMPS software package.

LAMMPS stands for Large-scale Atomic/Molecular Massively Parallel
Simulator.

Copyright (2003) Sandia Corporation.  Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software.  This software is distributed under
the GNU General Public License.

----------------------------------------------------------------------

LAMMPS is a classical molecular dynamics simulation code designed to
run efficiently on parallel computers.  It was developed at Sandia
National Laboratories, a US Department of Energy facility, with
funding from the DOE.  It is an open-source code, distributed freely
under the terms of the GNU Public License (GPL) version 2.

The code is maintained by the LAMMPS development team who can be emailed
at developers@lammps.org.  The LAMMPS WWW Site at www.lammps.org has
more information about the code and its uses.

The LAMMPS distribution includes the following files and directories:

README                     this file
LICENSE                    the GNU General Public License (GPLv2)
CITATION.cff               Citation information for LAMMPS in CFF format
bench                      benchmark inputs
cmake                      CMake build files
doc                        documentation
examples                   example inputs for many LAMMPS commands
fortran                    Fortran 2003 module for LAMMPS
lib                        additional provided or external libraries
potentials                 interatomic potential files
python                     Python module for LAMMPS
src                        source files
third_party                Copies of thirdparty software bundled with LAMMPS
tools                      pre- and post-processing tools
unittest                   test programs for use with CTest
.github                    Git and GitHub related files and tools

Point your browser at any of these files to get started:

https://docs.lammps.org/Manual.html         LAMMPS manual
https://docs.lammps.org/Intro.html          hi-level introduction
https://docs.lammps.org/Build.html          how to build LAMMPS
https://docs.lammps.org/Run_head.html       how to run LAMMPS
https://docs.lammps.org/Commands_all.html   Table of available commands
https://docs.lammps.org/Howto.html          Short tutorials and HowTo discussions
https://docs.lammps.org/Errors.html         How to interpret and debug errors
https://docs.lammps.org/Library.html        LAMMPS library interfaces
https://docs.lammps.org/Modify.html         how to modify and extend LAMMPS
https://docs.lammps.org/Developer.html      LAMMPS developer info

You can also create these doc pages locally:

% cd doc
% make html                # creates HTML pages in doc/html
% make pdf                 # creates Manual.pdf

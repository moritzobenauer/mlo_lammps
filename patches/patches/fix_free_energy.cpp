#include "fix_free_energy.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "group.h"
#include "modify.h"
#include "utils.h"
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

FixFreeEnergy::FixFreeEnergy(LAMMPS *_lmp, int narg, char **arg) : Fix(_lmp, narg, arg)
{
  // Expected syntax: fix ID group free_energy a b f [disable_reactions]
  if (narg < 6 || narg > 7)
    error->all(FLERR,
               "Illegal fix free_energy command. Expected: fix ID group free_energy a b f [disable_reactions]");

  scalar_flag = 1;
  extscalar = 1;
  energy_global_flag = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  // NOTE: a multiplies the quartic term only, it is not an overall scale of phi(z).
  //       phi(z) = a*z^4 - b*z^2 + f*z
  coeff_a = utils::numeric(FLERR, arg[3], false, lmp); // z^4 coeff
  coeff_b = utils::numeric(FLERR, arg[4], false, lmp); // z^2 coeff
  coeff_f = utils::numeric(FLERR, arg[5], false, lmp); // z coeff

  disable_reactions = 0;
  if (narg == 7) {
    if (strcmp(arg[6], "disable_reactions") == 0) {
      disable_reactions = 1;
    } else {
      error->all(FLERR,
                 "Illegal fix free_energy optional argument. Supported optional keyword: disable_reactions");
    }
  }

  // FIXED by CLAUDE (2026-08-31): with disable_reactions the z degree of freedom is
  // frozen, so advertise the removed DOF to temperature computes, which were reading
  // 2/3 of the true in-plane temperature. NOTE: this makes plain `compute temp` correct,
  // but `compute temp/partial 1 1 0` will then over-report by 3/2, because it prorates
  // fix-removed DOF evenly across its active dimensions (right for SHAKE-like isotropic
  // constraints, wrong for a purely-z freeze). Use `compute temp` in that case.
  if (disable_reactions) dof_flag = 1;

  e_total_all = 0.0;
  eflag = 0;
}

FixFreeEnergy::~FixFreeEnergy() {}

int FixFreeEnergy::setmask()
{
  int mask = 0;
  mask |= FixConst::POST_FORCE;
  mask |= FixConst::POST_FORCE_RESPA;
  // FIXED by CLAUDE (2026-08-31): phi(z) was absent during minimization, so the
  // pair style's z-force had nothing to oppose it and minimize drove all atoms to
  // lambda = 1.
  mask |= FixConst::MIN_POST_FORCE;
  // FIXED by CLAUDE (2026-08-31): keep v_z pinned outside of post_force (see below).
  if (disable_reactions) mask |= FixConst::POST_INTEGRATE;
  return mask;
}

void FixFreeEnergy::init()
{
  if (utils::strmatch(update->integrate_style, "^respa")) {
    ilevel_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels - 1;
  }

  // FIXED by CLAUDE (2026-08-31): disable_reactions only suppresses the z-forces of
  // other fixes (e.g. the Langevin z-kick) if this fix's post_force runs last.
  // That ordering used to be an undocumented, silent requirement - enforce it here.
  if (disable_reactions) {
    const auto &fixes = modify->get_fix_list();
    bool seen_self = false;
    for (const auto &ifix : fixes) {
      if (ifix == this) { seen_self = true; continue; }
      if (seen_self && (ifix->setmask() & FixConst::POST_FORCE))
        error->all(FLERR,
                   "fix free_energy disable_reactions requires that this fix be defined after "
                   "all other fixes that apply forces (offending later fix: {})", ifix->id);
    }
  }
}

void FixFreeEnergy::setup(int vflag)
{
  // FIXED by CLAUDE (2026-08-31): freeze v_z once here instead of every post_force
  // call - post_force is also invoked for energy-only re-evaluations (e.g. every
  // fix widom insertion trial), where mutating velocities corrupts the trajectory.
  if (disable_reactions) freeze_z_velocities();
  post_force(vflag);
}

void FixFreeEnergy::min_setup(int vflag) { setup(vflag); }

void FixFreeEnergy::freeze_z_velocities()
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) v[i][2] = 0.0;
}

void FixFreeEnergy::post_integrate() { freeze_z_velocities(); }

void FixFreeEnergy::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  eflag = 0;
  e_total = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double z = x[i][2];
      double z2 = z * z;

      if (disable_reactions) {
        // z is frozen: drop every z-force, including those added by earlier fixes.
        f[i][2] = 0.0;
      } else {
        // Force Fz = -dphi/dz = -(4*a*z^3 - 2*b*z + f)
        f[i][2] += -(4.0 * coeff_a * z2 * z - 2.0 * coeff_b * z + coeff_f);
      }

      // Energy phi(z) = a*z^4 - b*z^2 + f*z
      e_total += coeff_a * z2 * z2 - coeff_b * z2 + coeff_f * z;
    }
  }
}

void FixFreeEnergy::min_post_force(int vflag) { post_force(vflag); }

void FixFreeEnergy::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

// FIXED by CLAUDE (2026-08-31): number of DOF removed by freezing z, reported to
// temperature computes so that `compute temp` is not off by a factor of 3/2.
bigint FixFreeEnergy::dof(int igroup)
{
  if (!disable_reactions) return 0;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int other_bit = group->bitmask[igroup];

  bigint n = 0;
  for (int i = 0; i < nlocal; i++)
    if ((mask[i] & groupbit) && (mask[i] & other_bit)) n++;

  bigint nall = 0;
  MPI_Allreduce(&n, &nall, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  return nall;
}

double FixFreeEnergy::compute_scalar()
{
  // FIXED by CLAUDE (2026-08-31): the cache flag was reset on every call, so each
  // call did a fresh MPI_Allreduce. eflag is cleared in post_force instead.
  if (eflag == 0) {
    MPI_Allreduce(&e_total, &e_total_all, 1, MPI_DOUBLE, MPI_SUM, world);
    eflag = 1;
  }
  return e_total_all;
}

/* ----------------------------------------------------------------------
   expose the landscape coefficients of phi(z) = a*z^4 - b*z^2 + f*z

   pair_style lj/cut/mlo uses these to solve phi'(z) = 4a z^3 - 2b z + f = 0 for the
   barrier top, which is where Theta(z) must cross 0.5. Scalars, so dim = 0.
------------------------------------------------------------------------- */

void *FixFreeEnergy::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str, "a") == 0) return (void *) &coeff_a;
  if (strcmp(str, "b") == 0) return (void *) &coeff_b;
  if (strcmp(str, "f") == 0) return (void *) &coeff_f;
  return nullptr;
}

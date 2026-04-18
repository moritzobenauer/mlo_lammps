#include "fix_free_energy.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
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
  
  e_total_all = 0.0;
  eflag = 0;
}

FixFreeEnergy::~FixFreeEnergy() {}

int FixFreeEnergy::setmask()
{
  int mask = 0;
  mask |= FixConst::POST_FORCE;
  mask |= FixConst::POST_FORCE_RESPA;
  return mask;
}

void FixFreeEnergy::init()
{
  if (utils::strmatch(update->integrate_style, "^respa")) {
    ilevel_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels - 1;
  }
}

void FixFreeEnergy::setup(int vflag) { post_force(vflag); }

void FixFreeEnergy::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double **f = atom->f;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  eflag = 0;
  e_total = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double z = x[i][2];
      
      // Force Fz = -dV/dz = -(4*a*z^3 - 2*b*z + f)
      // Based on your prompt: V = z^4 - b*z^2 + f*z (assuming a=1)
      double f_z = -(4.0 * coeff_a * pow(z, 3) - 2.0 * coeff_b * z + coeff_f);
      f[i][2] += f_z;

      // Energy V = a*z^4 - b*z^2 + f*z
      e_total += coeff_a * pow(z, 4) - coeff_b * pow(z, 2) + coeff_f * z;

      if (disable_reactions) {
        v[i][2] = 0.0;
        f[i][2] = 0.0;
      }
    }
  }




}

void FixFreeEnergy::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

double FixFreeEnergy::compute_scalar()
{
  if (eflag == 0) {
    MPI_Allreduce(&e_total, &e_total_all, 1, MPI_DOUBLE, MPI_SUM, world);
    eflag = 1;
  }
  eflag = 0;
  return e_total_all;
}
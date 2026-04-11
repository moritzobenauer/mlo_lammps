/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Moritz L. Obenauer @ Princeton University, 2026
   lj/cut potential is only active in the xy plane, no z component of the force/energy
   Both force and energy are multiplied by a correction factor h(z) which is 1 for z>=0 and 0 for z<0

------------------------------------------------------------------------- */

#include "pair_lj_cut_mlo.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
// #include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

// To-Do: Disable respa

/* ---------------------------------------------------------------------- */

MLOPairLJCut::MLOPairLJCut(LAMMPS *lmp) : Pair(lmp)
{
  born_matrix_enable = 1;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

MLOPairLJCut::~MLOPairLJCut()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);
    memory->destroy(LJ_MINIMUM);
  }
}

/* ---------------------------------------------------------------------- */

void MLOPairLJCut::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq, r2inv, r6inv, forcelj, factor_lj, rsq_2d;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // Check if the current atom is active or inactive and calculate Theta(z)


    double sigmoid_alpha = 100.0;

    double sigmoid_arg = sigmoid_alpha * ztmp;    // Steepness of the sigmoid function
    double theta_i = 1.0 / (1.0 + exp(-sigmoid_arg));    // Smoothly transitions from 0 to 1 around z=0

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      // printf("Atom %d interacting with Atom %d: delx=%.3f, dely=%.3f, delz=%.3f, rsq=%.3f\n", i, j, delx, dely, delz, rsq);

      rsq_2d = delx * delx + dely * dely;    // Only consider xy components for cutoff

      jtype = type[j];

      double force_z = 0.0;


      // if (rsq < cutsq[itype][jtype]) {
      if (rsq_2d < cutsq[itype][jtype]) {



      // Calculate l based on the z position

        double sigmoid_arg = sigmoid_alpha * ztmp;    // Steepness of the sigmoid function
        double theta_j =
            1.0 / (1.0 + exp(-sigmoid_arg));    // Smoothly transitions from 0 to 1 around z=0

        // The force correction is Lambda = Theta(z_i) * Theta(z_j)
        double force_correction = theta_i * theta_j;
        
        // Debugging setting lambda = 1.0
        // force_correction = 0.0;
        // 

        // Calculating r^-1 and r^-6 for the 2D distances
        double r2inv_2d = 1.0 / rsq_2d;
        double r6inv_2d = r2inv_2d * r2inv_2d * r2inv_2d;


        // calculate the derivative for the potential in z \theta_i * d/dz \theta_j

        double force_derivative_factor_z = theta_i * (sigmoid_alpha * exp(-sigmoid_arg) / pow(1.0 + exp(-sigmoid_arg), 2));



        if (rsq_2d <= LJ_MINIMUM[itype][jtype]) {
          // AH force: return d/dr LJ
          forcelj = r6inv_2d * (lj1[itype][jtype] * r6inv_2d - lj2[itype][jtype]);

          force_z = epsilon[itype][jtype] * theta_i * force_derivative_factor_z;

        }

        else {
          // AH force: return l * d/dr LJ
          forcelj =
              r6inv_2d * (lj1[itype][jtype] * r6inv_2d - lj2[itype][jtype]) * force_correction;

          force_z = -((lj3[itype][jtype] * r6inv_2d - lj4[itype][jtype]) -
                offset[itype][jtype]) * theta_i * force_derivative_factor_z;
        }

        // printf("particle position: %f, force_z: %f, derivative: %f\n", x[i][2], force_z, force_derivative_factor_z);


        // Converting units and converting potential to force
        fpair = factor_lj * forcelj * r2inv_2d;

        double fpair_z = factor_lj * force_z;

        // Don't apply any forces in the z direction
        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += fpair_z;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= fpair_z;
        }

        if (eflag) {
         

          if (rsq_2d <= LJ_MINIMUM[itype][jtype]) {
            // AH: return LJ + (1-l) * epsilon
            evdwl = r6inv_2d * (lj3[itype][jtype] * r6inv_2d - lj4[itype][jtype]) +
                (1.0 - force_correction) * epsilon[itype][jtype];
          }
          else {
            // AH: return l * LJ
            evdwl = force_correction * (lj3[itype][jtype] * r6inv_2d - lj4[itype][jtype]) -
                offset[itype][jtype];
          }

          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
};

// Removed the entire respa stuff

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void MLOPairLJCut::allocate()
{
  allocated = 1;
  int n = atom->ntypes + 1;

  memory->create(setflag, n, n, "pair:setflag");
  for (int i = 1; i < n; i++)
    for (int j = i; j < n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n, n, "pair:cutsq");

  memory->create(cut, n, n, "pair:cut");
  memory->create(epsilon, n, n, "pair:epsilon");
  memory->create(sigma, n, n, "pair:sigma");
  memory->create(lj1, n, n, "pair:lj1");
  memory->create(lj2, n, n, "pair:lj2");
  memory->create(lj3, n, n, "pair:lj3");
  memory->create(lj4, n, n, "pair:lj4");
  memory->create(offset, n, n, "pair:offset");
  memory->create(LJ_MINIMUM, n, n, "pair:LJ_MINIMUM");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void MLOPairLJCut::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR, "Illegal pair_style command");

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void MLOPairLJCut::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5)
    error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double epsilon_one = utils::numeric(FLERR, arg[2], false, lmp);
  double sigma_one = utils::numeric(FLERR, arg[3], false, lmp);

  double cut_one = cut_global;
  if (narg == 5) cut_one = utils::numeric(FLERR, arg[4], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double MLOPairLJCut::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i], epsilon[j][j], sigma[i][i], sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i], sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i], cut[j][j]);
  }

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);

  double CUBE_ROOT_TWO = pow(2.0, 1.0/3.0);
  LJ_MINIMUM[i][j] = CUBE_ROOT_TWO * sigma[i][j] *sigma[i][j];    // Square Minimum of the LJ potential in 2D

  if (offset_flag && (cut[i][j] > 0.0)) {
    double ratio = sigma[i][j] / cut[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio, 12.0) - pow(ratio, 6.0));
  } else
    offset[i][j] = 0.0;

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  // check interior rRESPA cutoff

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2], all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count, all, 2, MPI_DOUBLE, MPI_SUM, world);

    double sig2 = sigma[i][j] * sigma[i][j];
    double sig6 = sig2 * sig2 * sig2;
    double rc3 = cut[i][j] * cut[i][j] * cut[i][j];
    double rc6 = rc3 * rc3;
    double rc9 = rc3 * rc6;
    double prefactor = 8.0 * MY_PI * all[0] * all[1] * epsilon[i][j] * sig6 / (9.0 * rc9);
    etail_ij = prefactor * (sig6 - 3.0 * rc6);
    ptail_ij = 2.0 * prefactor * (2.0 * sig6 - 3.0 * rc6);
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void MLOPairLJCut::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j], sizeof(double), 1, fp);
        fwrite(&sigma[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void MLOPairLJCut::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &epsilon[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &sigma[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&epsilon[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&sigma[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void MLOPairLJCut::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
  fwrite(&tail_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void MLOPairLJCut::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tail_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&tail_flag, 1, MPI_INT, 0, world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void MLOPairLJCut::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) fprintf(fp, "%d %g %g\n", i, epsilon[i][i], sigma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void MLOPairLJCut::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g\n", i, j, epsilon[i][j], sigma[i][j], cut[i][j]);
}

/* ---------------------------------------------------------------------- */

// double MLOPairLJCut::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
//                             double /*factor_coul*/, double factor_lj, double &fforce)
// {
//   double r2inv, r6inv, forcelj, philj;

//   r2inv = 1.0 / rsq;
//   r6inv = r2inv * r2inv * r2inv;
//   forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
//   fforce = factor_lj * forcelj * r2inv;

//   philj = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) - offset[itype][jtype];
//   return factor_lj * philj;
// }

/* ---------------------------------------------------------------------- */

void MLOPairLJCut::born_matrix(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                               double /*factor_coul*/, double factor_lj, double &dupair,
                               double &du2pair)
{
  double rinv, r2inv, r6inv, du, du2;

  r2inv = 1.0 / rsq;
  rinv = sqrt(r2inv);
  r6inv = r2inv * r2inv * r2inv;

  // Reminder: lj1 = 48*e*s^12, lj2 = 24*e*s^6
  // so dupair = -forcelj/r = -fforce*r (forcelj from single method)

  du = r6inv * rinv * (lj2[itype][jtype] - lj1[itype][jtype] * r6inv);
  du2 = r6inv * r2inv * (13 * lj1[itype][jtype] * r6inv - 7 * lj2[itype][jtype]);

  dupair = factor_lj * du;
  du2pair = factor_lj * du2;
}

/* ---------------------------------------------------------------------- */

void *MLOPairLJCut::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "epsilon") == 0) return (void *) epsilon;
  if (strcmp(str, "sigma") == 0) return (void *) sigma;
  return nullptr;
}

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
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
// #include "respa.h"
#include "update.h"

#include <cmath>
#include <cstdint>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

// "MLO_PAIR" -- leading magic that makes a pre-versioning restart file detectable.
static constexpr uint64_t MLO_RESTART_MAGIC = 0x4D4C4F5F50414952ULL;

// To-Do: Disable respa

/* ---------------------------------------------------------------------- */

MLOPairLJCut::MLOPairLJCut(LAMMPS *lmp) : Pair(lmp)
{
  // FIXED by CLAUDE (2026-08-31): single() and born_matrix() do not implement this
  // potential (xy-only distance, Theta(z) coupling), so they must not be advertised.
  // A stub returning a plausible number is worse than no support at all: fix widom and
  // fix gcmc silently promote themselves to the full-energy path when single_enable == 0.
  single_enable = 0;
  born_matrix_enable = 0;
  writedata = 1;

  // Sentinels: these settings are not stored in restart files, see init_style().
  cut_global = 0.0;
  z_star = 0.0;
  cutoff_2d = -1.0;
  sigmoid_alpha = 50.0;
  z_budget = -1.0;
  zbudget_explicit = 0;
  settings_set = 0;
  lambda_mix_flag = 0;
  zstar_auto = 1;
  zcheck_flag = 1;
  zperiodic_ok = 0;
  theta_check_on = 1;
  theta_tol = 1.0e-4;
  last_zcheck = -1;
  zwarn_issued = 0;
  restart_version = 2;
  z_star_type = nullptr;
}

/* ---------------------------------------------------------------------- */

MLOPairLJCut::~MLOPairLJCut()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(cut_2d);
    memory->destroy(cut_2d_sq);
    memory->destroy(lambda_max);
    memory->destroy(z_star_type);
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
  double xtmp, ytmp, delx, dely, evdwl, fpair;
  double forcelj, factor_lj, rsq_2d;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  check_z_budget();

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

  // FIXED by CLAUDE (2026-08-31): the in-plane cutoff is now actually applied. The gate
  // used to test the 3d cutsq (= cut_global^2) against a distance whose z component had
  // been zeroed, so the effective xy interaction range was cut_global (5.59 sigma in the
  // input scripts) instead of the declared cutoff_2d (2.5 sigma), and cutoff_2d was a
  // dead parameter. cut[i][j] remains the *neighbor* cutoff, which stays inflated so that
  // no in-plane neighbor is culled by its z separation (derived in init_one from z_budget).
  //
  // The gate is now per type pair: cut_2d_sq[itype][jtype], not one global scalar. Note
  // that cutsq[i][j] is NOT usable here -- Pair::init() builds it from init_one()'s return
  // value, which has to be the inflated neighbor cutoff.

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // Theta(z_i) is the smooth indicator of the bonding state. It has to be a functional
    // of the free energy landscape, because the position of the barrier top shifts
    // slightly with the landscape parameters.
    // Clamping the argument avoids overflow in either limit; Theta' is ~1e-25 at the
    // clamp, so nothing physical is discarded.

    double sigmoid_arg = sigmoid_alpha * (x[i][2] - z_star_type[itype]);
    sigmoid_arg = std::max(-60.0, std::min(60.0, sigmoid_arg));
    double exp_i = exp(-sigmoid_arg);
    double theta_i = 1.0 / (1.0 + exp_i);
    // Theta'(z) = alpha * e^-s / (1 + e^-s)^2
    double dtheta_i = sigmoid_alpha * exp_i * theta_i * theta_i;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      // The LJ interaction lives in the xy plane only: there is no z component of the
      // pair separation. z enters exclusively through Theta(z_i) * Theta(z_j).

      rsq_2d = delx * delx + dely * dely;
      jtype = type[j];

      if (rsq_2d >= cut_2d_sq[itype][jtype]) continue;

      // FIXED by CLAUDE (2026-08-31): guard the reciprocal (it used to be evaluated for
      // every neighbor before any cutoff test, so a coincident pair produced inf/NaN
      // forces silently). Zero in-plane distance is unphysical here - it also happens if
      // a particle sees its own periodic z image, see init_style().
      if (rsq_2d == 0.0)
        error->one(FLERR, "Pair style lj/cut/mlo: zero in-plane distance between atoms");

      // Exact early-outs. Both skip pairs whose energy, in-plane force and z-force are all
      // identically zero, so they change nothing numerically -- they only avoid the exp()
      // for Theta(z_j) below.
      //   epsilon == 0 : lj1..lj4 and offset are all zero, so both branches vanish.
      //   lambda_max == 0 outside r_m : the attractive branch is lambda*(...) with
      //     lambda = 0, and dU/dlambda is multiplied by lambda_max in the z-force.
      // Inside r_m a lambda_max == 0 pair still has its full WCA core, so it must NOT be
      // skipped there.
      if (epsilon[itype][jtype] == 0.0) continue;
      if (lambda_max[itype][jtype] == 0.0 && rsq_2d > LJ_MINIMUM[itype][jtype]) continue;

      double r2inv_2d = 1.0 / rsq_2d;
      double r6inv_2d = r2inv_2d * r2inv_2d * r2inv_2d;

      double sigmoid_arg_j = sigmoid_alpha * (x[j][2] - z_star_type[jtype]);
      sigmoid_arg_j = std::max(-60.0, std::min(60.0, sigmoid_arg_j));
      double exp_j = exp(-sigmoid_arg_j);
      double theta_j = 1.0 / (1.0 + exp_j);
      double dtheta_j = sigmoid_alpha * exp_j * theta_j * theta_j;

      // lambda_ij = lambda_max[i][j] * Theta(z_i) * Theta(z_j).
      //
      // The PRODUCT of the two indicators is deliberate, and is not the arithmetic mixing
      // rule of the published Ashbaugh-Hatch model: the pair attracts only if *both*
      // partners are in the bonding state, so a single inactive partner removes the
      // attraction for the pair. That is what makes "inactive states never attract" a
      // structural property rather than a parameter choice.
      //
      // lambda_max is the per-type-pair attraction scale, kept separate from epsilon so
      // that a pair can keep its full excluded volume with no attraction at all.
      double lam_max = lambda_max[itype][jtype];
      double lambda = lam_max * theta_i * theta_j;

      double u_lj = r6inv_2d * (lj3[itype][jtype] * r6inv_2d - lj4[itype][jtype]);
      double shift = offset[itype][jtype];

      // Ashbaugh-Hatch form, with the energy shift applied consistently in both branches
      // (FIXED by CLAUDE, 2026-08-31 - the repulsive branch carried no shift at all, and
      // the attractive branch shifted the energy by `offset` but its force by
      // `lambda*offset`, so the force was not the gradient of the energy):
      //
      //   r <= 2^(1/6) sigma : U = U_LJ + (1 - lambda) * eps - lambda * shift
      //   r >  2^(1/6) sigma : U = lambda * (U_LJ - shift)
      //
      // Both are continuous at the crossover and vanish at r = cutoff_2d, and the z-force
      // is -(dU/dlambda) * dlambda/dz with dlambda/dz_i = Theta_j * Theta_i'.

      double dU_dlambda;

      if (rsq_2d <= LJ_MINIMUM[itype][jtype]) {
        // AH force: full d/dr LJ
        forcelj = r6inv_2d * (lj1[itype][jtype] * r6inv_2d - lj2[itype][jtype]);
        dU_dlambda = -(epsilon[itype][jtype] + shift);
        if (eflag)
          evdwl = u_lj + (1.0 - lambda) * epsilon[itype][jtype] - lambda * shift;
      } else {
        // AH force: lambda * d/dr LJ
        forcelj = r6inv_2d * (lj1[itype][jtype] * r6inv_2d - lj2[itype][jtype]) * lambda;
        dU_dlambda = u_lj - shift;
        if (eflag) evdwl = lambda * (u_lj - shift);
      }

      // d(lambda)/dz_i = lambda_max * Theta_j * Theta_i', so lambda_max carries through to
      // the reaction-coordinate forces as well. Dropping it here while keeping it in
      // `lambda` above would leave the force no longer the gradient of the energy.
      double force_z_i = -dU_dlambda * lam_max * theta_j * dtheta_i;
      double force_z_j = -dU_dlambda * lam_max * theta_i * dtheta_j;

      // Converting units and converting potential to force
      fpair = factor_lj * forcelj * r2inv_2d;

      f[i][0] += delx * fpair;
      f[i][1] += dely * fpair;
      f[i][2] += factor_lj * force_z_i;

      if (newton_pair || j < nlocal) {
        // x and y are cyclic coordinates --> momentum conservation requires f_j = -f_i
        f[j][0] -= delx * fpair;
        f[j][1] -= dely * fpair;
        // z is not a cyclic coordinate --> no momentum conservation. U depends on z_i and
        // z_j separately rather than on z_i - z_j, so f_j gets its own Theta_j' term.
        f[j][2] += factor_lj * force_z_j;
      }

      if (eflag) evdwl *= factor_lj;

      // NOTE: delz = 0 is passed on purpose, so ev_tally builds the correct xy pair
      // virial. The global virial nevertheless comes from virial_fdotr_compute() below,
      // which sums x.f over the force array and therefore also picks up the
      // non-pairwise Theta(z) z-forces. On top of that the box volume contains Lz.
      // Do not quote `thermo press` for this pair style; compute a 2d virial explicitly
      // if pressure is needed.
      if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, 0.0);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

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
  memory->create(cut_2d, n, n, "pair:cut_2d");
  memory->create(cut_2d_sq, n, n, "pair:cut_2d_sq");
  memory->create(lambda_max, n, n, "pair:lambda_max");
  memory->create(z_star_type, n, "pair:z_star_type");
  for (int i = 0; i < n; i++) z_star_type[i] = z_star;
  memory->create(epsilon, n, n, "pair:epsilon");
  memory->create(sigma, n, n, "pair:sigma");
  memory->create(lj1, n, n, "pair:lj1");
  memory->create(lj2, n, n, "pair:lj2");
  memory->create(lj3, n, n, "pair:lj3");
  memory->create(lj4, n, n, "pair:lj4");
  memory->create(offset, n, n, "pair:offset");
  memory->create(LJ_MINIMUM, n, n, "pair:LJ_MINIMUM");

  // Initialize the FULL square, not just the upper triangle: compute() indexes
  // [itype][jtype] with no ordering guarantee, and read_restart of a pre-cut_2d file
  // relies on these defaults. -1.0 means "not set yet", resolved in init_one().
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      cut_2d[i][j] = -1.0;
      cut_2d_sq[i][j] = 0.0;
      lambda_max[i][j] = 1.0;    // full attraction: reproduces the pre-lambda_max potential
    }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void MLOPairLJCut::settings(int narg, char **arg)
{
  // pair_style lj/cut/mlo <neighbor cutoff> <z_star> <xy cutoff> [alpha] [keyword values...]
  if (narg < 3)
    error->all(FLERR,
               "Illegal pair_style lj/cut/mlo command. Expected: pair_style lj/cut/mlo "
               "cutoff_3d z_star cutoff_xy [alpha] [zbudget dz]");

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);
  z_star = utils::numeric(FLERR, arg[1], false, lmp);
  cutoff_2d = utils::numeric(FLERR, arg[2], false, lmp);

  zbudget_explicit = 0;

  // FIXED by CLAUDE (2026-08-31): the steepness of Theta(z) was hardcoded (and had drifted
  // to a different value in the single() code path), while z_star was already an argument.
  // Optional, so existing 3-argument input scripts keep the previous alpha = 50.
  //
  // The legacy 4th positional argument is still accepted, so keyword parsing starts at
  // arg[3] only when arg[3] does not look like a number.
  int iarg = 3;
  if (narg > 3 && utils::is_double(arg[3])) {
    sigmoid_alpha = utils::numeric(FLERR, arg[3], false, lmp);
    iarg = 4;
  }

  while (iarg < narg) {
    if (strcmp(arg[iarg], "alpha") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "pair_style lj/cut/mlo alpha", error);
      sigmoid_alpha = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "zbudget") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "pair_style lj/cut/mlo zbudget", error);
      z_budget = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      zbudget_explicit = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg], "zstar") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "pair_style lj/cut/mlo zstar", error);
      if (strcmp(arg[iarg + 1], "auto") == 0)
        zstar_auto = 1;
      else if (strcmp(arg[iarg + 1], "fixed") == 0)
        zstar_auto = 0;
      else
        error->all(FLERR, "Pair style lj/cut/mlo zstar must be auto or fixed");
      iarg += 2;
    } else if (strcmp(arg[iarg], "zcheck") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "pair_style lj/cut/mlo zcheck", error);
      if (strcmp(arg[iarg + 1], "off") == 0)
        zcheck_flag = 0;
      else if (strcmp(arg[iarg + 1], "warn") == 0)
        zcheck_flag = 1;
      else if (strcmp(arg[iarg + 1], "strict") == 0)
        zcheck_flag = 2;
      else
        error->all(FLERR, "Pair style lj/cut/mlo zcheck must be off, warn or strict");
      iarg += 2;
    } else if (strcmp(arg[iarg], "theta_check") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, "pair_style lj/cut/mlo theta_check", error);
      if (strcmp(arg[iarg + 1], "off") == 0) {
        theta_check_on = 0;
      } else {
        theta_tol = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
        theta_check_on = 1;
        if (theta_tol <= 0.0)
          error->all(FLERR, "Pair style lj/cut/mlo theta_check tolerance must be > 0.0");
      }
      iarg += 2;
    } else if (strcmp(arg[iarg], "zperiodic_ok") == 0) {
      zperiodic_ok = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg], "lambda_mix") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, "pair_style lj/cut/mlo lambda_mix", error);
      if (strcmp(arg[iarg + 1], "geometric") == 0)
        lambda_mix_flag = 0;
      else if (strcmp(arg[iarg + 1], "arithmetic") == 0)
        lambda_mix_flag = 1;
      else
        error->all(FLERR, "Pair style lj/cut/mlo lambda_mix must be geometric or arithmetic");
      iarg += 2;
    } else {
      error->all(FLERR, "Unknown pair_style lj/cut/mlo keyword: {}", arg[iarg]);
    }
  }

  if (cutoff_2d <= 0.0) error->all(FLERR, "Pair style lj/cut/mlo xy cutoff must be > 0.0");
  if (sigmoid_alpha <= 0.0) error->all(FLERR, "Pair style lj/cut/mlo alpha must be > 0.0");
  if (cutoff_2d > cut_global)
    error->all(FLERR,
               "Pair style lj/cut/mlo xy cutoff ({}) must not exceed the neighbor cutoff "
               "({}), otherwise interacting pairs are missing from the neighbor list",
               cutoff_2d, cut_global);

  // The three cutoff-ish quantities have only two degrees of freedom, because
  //   cut_neigh = sqrt(cut_2d^2 + z_budget^2).
  // The legacy signature exposes the two derived ones and leaves z_budget implicit, which
  // is exactly the `cutoff = sqrt(cutoff_xy^2 + z_max^2)` arithmetic the input scripts do
  // by hand. Recover z_budget from them so old scripts are unchanged, and let `zbudget`
  // override it.
  if (!zbudget_explicit) z_budget = sqrt(cut_global * cut_global - cutoff_2d * cutoff_2d);

  if (z_budget <= 0.0)
    error->all(FLERR,
               "Pair style lj/cut/mlo has a zero z budget (neighbor cutoff {} equals the xy "
               "cutoff {}), so any z separation culls an in-plane neighbor and the pair "
               "silently loses even its repulsive core. Increase the neighbor cutoff or set "
               "the zbudget keyword",
               cut_global, cutoff_2d);

  settings_set = 1;

  if (allocated)
    for (int i = 0; i <= atom->ntypes; i++) z_star_type[i] = z_star;

  // reset in-plane cutoffs that have been explicitly set

  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut_2d[i][j] = cutoff_2d;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void MLOPairLJCut::coeff(int narg, char **arg)
{
  if (narg < 4) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double epsilon_one = utils::numeric(FLERR, arg[2], false, lmp);
  double sigma_one = utils::numeric(FLERR, arg[3], false, lmp);

  // The in-plane cutoff is now per type pair, and it is a KEYWORD, not a positional.
  // A bare 5th positional used to set cut[i][j], i.e. the deliberately inflated *neighbor*
  // cutoff, not the interaction range -- and the old write_data_all() emitted exactly that
  // value in that slot. Silently reinterpreting it as an in-plane cutoff would reproduce
  // the 5.59-sigma interaction range of CODE_REVIEW P1 #1, so reject it loudly instead.
  double cut_2d_one = (cutoff_2d > 0.0) ? cutoff_2d : -1.0;
  double lambda_one = 1.0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "lambda") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "pair_coeff lj/cut/mlo lambda", error);
      lambda_one = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "cut") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "pair_coeff lj/cut/mlo cut", error);
      cut_2d_one = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (utils::is_double(arg[iarg])) {
      error->all(FLERR,
                 "pair_coeff for lj/cut/mlo no longer takes a positional 5th argument. It "
                 "used to set the (inflated 3D) neighbor cutoff, which is now derived from "
                 "the zbudget keyword. Write `cut {}` to set the in-plane interaction "
                 "cutoff explicitly",
                 arg[iarg]);
    } else {
      error->all(FLERR, "Unknown pair_coeff lj/cut/mlo keyword: {}", arg[iarg]);
    }
  }

  if (epsilon_one < 0.0) error->all(FLERR, "Pair style lj/cut/mlo epsilon must be >= 0.0");
  if (sigma_one <= 0.0) error->all(FLERR, "Pair style lj/cut/mlo sigma must be > 0.0");
  if (cut_2d_one <= 0.0)
    error->all(FLERR, "Pair style lj/cut/mlo in-plane cutoff must be > 0.0");
  if (lambda_one < 0.0) error->all(FLERR, "Pair style lj/cut/mlo lambda must be >= 0.0");
  // lambda > 1 is well posed -- it just means the well is deeper than the core amplitude,
  // still continuous at r_m and still a gradient field -- but the name implies [0,1].
  if (lambda_one > 1.0 && comm->me == 0)
    error->warning(FLERR,
                   "Pair style lj/cut/mlo lambda {} exceeds 1.0: the attractive well is "
                   "deeper than epsilon",
                   lambda_one);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut_2d[i][j] = cut_2d_one;
      lambda_max[i][j] = lambda_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
}

/* ----------------------------------------------------------------------
   init specific to this pair style

   FIXED by CLAUDE (2026-08-31): the geometric assumptions of this pair style used to be
   unstated and unchecked. Because the interaction is gated on the in-plane distance alone
   while the neighbor list is built from true 3d distances, a pair that is close in xy but
   far apart in z drops out of the neighbor list and loses even its repulsive core.
------------------------------------------------------------------------- */

void MLOPairLJCut::init_style()
{
  Pair::init_style();

  // These settings are not stored in restart files, so a pair_style command has to be
  // reissued after read_restart. Fail loudly instead of running with undefined values.
  if (!settings_set)
    error->all(FLERR,
               "Pair style lj/cut/mlo has no pair_style settings. They are not stored in "
               "restart files - reissue pair_style after read_restart");

  if (tail_flag)
    error->all(FLERR,
               "Pair style lj/cut/mlo does not support pair_modify tail yes: the 3d "
               "tail-correction formula does not apply to an xy-only interaction");

  // (a0) z must not be domain-decomposed. With procgrid[2] == 1 and non-periodic z,
  //      Comm::setup gives maxneed[2] == 0: there is no z-direction ghost exchange, ghosts
  //      are selected on xy proximity alone, and a rank's local+ghost set therefore holds
  //      every atom that could be an in-plane neighbor of one of its locals at ANY z. That
  //      is what makes the per-rebuild guard in check_z_budget() exact and collective-free.
  //      It also avoids a badly imbalanced decomposition, since every atom sits near z = 0.

  if (comm->procgrid[2] != 1)
    error->all(FLERR,
               "Pair style lj/cut/mlo requires a single processor along z (procgrid is {}), "
               "because z is a chemical reaction coordinate rather than a spatial dimension. "
               "Add `processors * * 1`",
               comm->procgrid[2]);

  // (a1) periodic z is hazardous, and not in a way the old Lz > 2*cutneigh check covered.
  //      AtomVec::pack_comm stores a periodic-image ghost at x[j][2] + pbc[2]*zprd, and
  //      compute() reads x[j][2] to evaluate Theta(z_j). So a z-image ghost of an atom in
  //      the inactive well is seen at z + Lz, where Theta ~ 1 instead of ~0: its chemical
  //      state, and with it the sign of the attraction, silently inverts. Wrapping z is
  //      just as bad -- Theta jumps discontinuously and phi(z) delivers a delta-function
  //      force. A non-periodic z removes the whole class.

  if (domain->zperiodic && !zperiodic_ok)
    error->all(FLERR,
               "Pair style lj/cut/mlo requires a non-periodic z (`boundary p p f`): z is a "
               "chemical reaction coordinate, and a periodic image ghost carries z + Lz, so "
               "Theta(z) evaluates its chemical state at the wrong z. Use `boundary p p f`, "
               "or add the zperiodic_ok keyword to accept the risk");

  // (a2) the z spread must actually fit the budget. The previous version of this check
  //     sampled the instantaneous z spread and treated it as the requirement; that is
  //     backwards -- z is a dynamical variable, so the budget is the contract and the
  //     spread is what must be held inside it. A pair close in xy but separated by more
  //     than z_budget in z is absent from the neighbor list and loses even its repulsive
  //     core, silently.
  //
  //     The complementary guarantee -- that each pair's neighbor cutoff is wide enough for
  //     the budget -- is asserted in init_one(), which is where cut[i][j] is derived.
  //     Pair::init() runs init_style() BEFORE the init_one() loop, so cut[i][j] and the
  //     mixed cut_2d[i][j] do not exist yet at this point.

  double **x = atom->x;
  int nlocal = atom->nlocal;
  double zlo = 0.0, zhi = 0.0;
  if (nlocal > 0) zlo = zhi = x[0][2];
  for (int i = 1; i < nlocal; i++) {
    zlo = MIN(zlo, x[i][2]);
    zhi = MAX(zhi, x[i][2]);
  }
  double zlo_all, zhi_all;
  MPI_Allreduce(&zlo, &zlo_all, 1, MPI_DOUBLE, MPI_MIN, world);
  MPI_Allreduce(&zhi, &zhi_all, 1, MPI_DOUBLE, MPI_MAX, world);

  double dz_max = zhi_all - zlo_all;
  if (domain->zperiodic) dz_max = MIN(dz_max, 0.5 * domain->zprd);
  if (dz_max > z_budget)
    error->all(FLERR,
               "Pair style lj/cut/mlo: the z spread is {} but the z budget is only {}. "
               "In-plane neighbors separated by more than the budget in z are culled from "
               "the neighbor list and lose even their repulsive core. Raise the neighbor "
               "cutoff (first pair_style argument) or set `zbudget {}`",
               dz_max, z_budget, dz_max);

  // Tie the Theta(z) crossover to the actual free-energy landscape, and report how close
  // to zero the residual inactive-inactive attraction really is.
  derive_z_star();
  if (theta_check_on) report_theta_residual();

  // (b) with periodic z a particle must never see its own z image, which sits at zero
  //     in-plane distance and would give an infinite force.

  double cutneigh = cut_global + neighbor->skin;
  if (domain->zperiodic && domain->zprd <= 2.0 * cutneigh)
    error->all(FLERR,
               "Pair style lj/cut/mlo requires Lz > 2*(cutoff+skin) = {} with periodic z, "
               "so that no particle interacts with its own z image at zero in-plane "
               "distance. Lz is currently {}",
               2.0 * cutneigh, domain->zprd);
}

/* ----------------------------------------------------------------------
   per-rebuild guard on the z spread

   The interaction is gated on the in-plane distance, but the neighbor list is built from
   true 3d distances. A pair that is close in xy yet separated by more than z_budget in z is
   therefore absent from the list and loses even its repulsive core - silently, with a
   plausible-looking trajectory. init_style() checks this once; z is dynamical, so it has to
   be rechecked whenever the list is rebuilt.

   Cost: with procgrid[2] == 1 and non-periodic z (both enforced in init_style), Comm::setup
   leaves maxneed[2] == 0, so ghosts are chosen on xy proximity alone and a rank's
   local+ghost set contains every possible in-plane neighbor of its locals at any z. The
   rank-local spread is then the exact bound for that rank's list, and NO collective is
   needed. The periodic fallback costs one 2-double Allreduce per rebuild.
------------------------------------------------------------------------- */

void MLOPairLJCut::check_z_budget()
{
  if (zcheck_flag == 0) return;
  if (neighbor->ago != 0) return;
  if (update->ntimestep == last_zcheck) return;
  last_zcheck = update->ntimestep;

  double **x = atom->x;
  int nall = atom->nlocal;
  if (!domain->zperiodic) nall += atom->nghost;
  if (nall == 0) return;

  double zlo = x[0][2], zhi = x[0][2];
  for (int k = 1; k < nall; k++) {
    zlo = MIN(zlo, x[k][2]);
    zhi = MAX(zhi, x[k][2]);
  }
  double dz = zhi - zlo;

  if (domain->zperiodic) {
    double in[2] = {zhi, -zlo}, out[2];
    MPI_Allreduce(in, out, 2, MPI_DOUBLE, MPI_MAX, world);
    dz = MIN(out[0] + out[1], 0.5 * domain->zprd);
  }

  if (dz <= z_budget) return;

  if (zcheck_flag == 2)
    error->one(FLERR,
               "Pair style lj/cut/mlo: z spread {} exceeds the z budget {} at step {}. "
               "In-plane neighbors separated by more than the budget in z are culled from "
               "the neighbor list and lose even their repulsive core",
               dz, z_budget, update->ntimestep);
  else if (!zwarn_issued) {
    zwarn_issued = 1;
    error->warning(FLERR,
                   "Pair style lj/cut/mlo: z spread {} exceeds the z budget {} at step {}; "
                   "in-plane neighbors separated by more than the budget in z lose even "
                   "their repulsive core. Raise zbudget or confine z. Warning issued once",
                   dz, z_budget, update->ntimestep);
  }
}

/* ----------------------------------------------------------------------
   extrema of phi(z) = a z^4 - b z^2 + f z

   Solves phi'(z) = 4a z^3 - 2b z + f = 0 by bisection rather than Cardano: phi' has its own
   extrema at z = +/- sqrt(b/6a), which bracket the three roots exactly, and bisection cannot
   lose precision to catastrophic cancellation the way the trigonometric form can near a
   double root. Runs once per run.

   Returns 1 with zlo < ztop < zhi (inactive minimum, barrier top, bonding minimum) when phi
   is a genuine double well, 0 otherwise.
------------------------------------------------------------------------- */

static double dphi(double a, double b, double f, double z)
{
  return 4.0 * a * z * z * z - 2.0 * b * z + f;
}

int MLOPairLJCut::landscape_extrema(double a, double b, double f, double &zlo, double &ztop,
                                    double &zhi)
{
  if (a <= 0.0 || b <= 0.0) return 0;

  const double s = sqrt(b / (6.0 * a));
  if (dphi(a, b, f, -s) <= 0.0 || dphi(a, b, f, s) >= 0.0) return 0;

  const double Z = 2.0 * sqrt(b / (2.0 * a)) + 1.0;
  double bracket[3][2] = {{-Z, -s}, {-s, s}, {s, Z}};
  double root[3];

  for (int k = 0; k < 3; k++) {
    double x0 = bracket[k][0], x1 = bracket[k][1];
    double f0 = dphi(a, b, f, x0);
    if (dphi(a, b, f, x0) * dphi(a, b, f, x1) > 0.0) return 0;
    for (int it = 0; it < 200; it++) {
      double xm = 0.5 * (x0 + x1);
      double fm = dphi(a, b, f, xm);
      if (fm == 0.0 || 0.5 * (x1 - x0) < 1.0e-15 * (1.0 + fabs(xm))) {
        x0 = x1 = xm;
        break;
      }
      if ((f0 < 0.0) == (fm < 0.0)) {
        x0 = xm;
        f0 = fm;
      } else
        x1 = xm;
    }
    root[k] = 0.5 * (x0 + x1);
  }

  zlo = root[0];
  ztop = root[1];
  zhi = root[2];
  return 1;
}

/* ----------------------------------------------------------------------
   map each atom type to its fix free_energy landscape, and put the Theta(z) crossover at
   that landscape's barrier top

   Theta has to be a functional of the free-energy landscape: the barrier top moves with b
   and f, and the landscape is already per type in practice (two_particles.in runs two fix
   free_energy instances on type-based groups, whose barrier tops are 0.132 and 0). A single
   hand-supplied z_star cannot sit on both.

   The type -> fix map is built from the atoms rather than by assuming the group was defined
   by type, so it stays correct for a group defined by region or by id.
------------------------------------------------------------------------- */

void MLOPairLJCut::derive_z_star()
{
  const int ntypes = atom->ntypes;
  for (int t = 0; t <= ntypes; t++) z_star_type[t] = z_star;

  auto fixes = modify->get_fix_by_style("free_energy");
  const int nfix = (int) fixes.size();
  if (nfix == 0) {
    if (comm->me == 0 && theta_check_on)
      utils::logmesg(lmp,
                     "lj/cut/mlo: no fix free_energy found; Theta(z) crosses at the supplied "
                     "z_star = {} for every type\n",
                     z_star);
    return;
  }
  if (nfix > 31)
    error->all(FLERR, "Pair style lj/cut/mlo supports at most 31 fix free_energy instances");

  // which landscapes claim each type (bitmask, OR-reduced across ranks)
  int *claim = new int[ntypes + 1]();
  int *mask = atom->mask;
  int *type = atom->type;
  for (int k = 0; k < nfix; k++)
    for (int i = 0; i < atom->nlocal; i++)
      if (mask[i] & fixes[k]->groupbit) claim[type[i]] |= (1 << k);

  int *claim_all = new int[ntypes + 1];
  MPI_Allreduce(claim, claim_all, ntypes + 1, MPI_INT, MPI_BOR, world);

  int uncovered = 0;
  for (int t = 1; t <= ntypes; t++) {
    int bits = claim_all[t];
    int n = 0, which = -1;
    for (int k = 0; k < nfix; k++)
      if (bits & (1 << k)) {
        n++;
        which = k;
      }

    if (n == 0) {
      uncovered++;
      continue;    // no atoms of this type, or no landscape: keep the fallback
    }
    if (n > 1)
      error->all(FLERR,
                 "Pair style lj/cut/mlo: atom type {} is claimed by {} different fix "
                 "free_energy instances, so its Theta(z) crossover is ambiguous",
                 t, n);

    int dim;
    double *pa = (double *) fixes[which]->extract("a", dim);
    double *pb = (double *) fixes[which]->extract("b", dim);
    double *pf = (double *) fixes[which]->extract("f", dim);
    if (!pa || !pb || !pf)
      error->all(FLERR,
                 "Pair style lj/cut/mlo could not read the landscape coefficients from fix "
                 "{}",
                 fixes[which]->id);

    double zl, zt, zh;
    if (!landscape_extrema(*pa, *pb, *pf, zl, zt, zh)) {
      if (comm->me == 0)
        error->warning(FLERR,
                       "Pair style lj/cut/mlo: fix {} (a={}, b={}, f={}) is not a double "
                       "well, so there is no two-state interpretation for type {}; keeping "
                       "z_star = {}",
                       fixes[which]->id, *pa, *pb, *pf, t, z_star);
      continue;
    }
    if (zstar_auto) z_star_type[t] = zt;
  }

  if (uncovered && comm->me == 0)
    error->warning(FLERR,
                   "Pair style lj/cut/mlo: {} atom type(s) are not covered by any fix "
                   "free_energy. Such atoms have no phi(z), and the pair z-force is >= 0 "
                   "identically, so their z is unbounded and will run out of the z budget",
                   uncovered);

  delete[] claim;
  delete[] claim_all;
}

/* ----------------------------------------------------------------------
   report how close to zero the residual inactive-inactive attraction actually is

   lambda_ij = lambda_max * Theta_i * Theta_j is never exactly zero, so "inactive states
   never attract" is a statement about magnitude. At alpha = 50 the residual is ~1e-54 - far
   below double precision - but alpha is a user-facing argument now, and at alpha = 2 the
   residual is ~6e-3 of epsilon per neighbor. Report it, and say so when it is not negligible.

   Deliberately NOT fixed by clamping Theta to zero: a hard floor makes Theta' discontinuous,
   so U is no longer C1, symplectic energy conservation degrades to O(1) drift per barrier
   crossing, and the Langevin stationary distribution is biased. Worse, the kink would sit
   exactly where particles cross it ~k times per unit time, so the error would scale with the
   reaction rate being measured. Use theta_form switch (a C2 compactly-supported switch) if an
   identically-zero residual is needed.
------------------------------------------------------------------------- */

void MLOPairLJCut::report_theta_residual()
{
  auto fixes = modify->get_fix_by_style("free_energy");
  if (fixes.empty()) return;

  double lam_max_all = 0.0;
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) lam_max_all = MAX(lam_max_all, lambda_max[i][j]);

  double worst_res = 0.0;
  double worst_disp = 0.0;

  if (comm->me == 0)
    utils::logmesg(lmp, "lj/cut/mlo state diagnostics (alpha = {}, crossover width 1/alpha = {}):\n",
                   sigmoid_alpha, 1.0 / sigmoid_alpha);

  for (auto &fx : fixes) {
    int dim;
    double *pa = (double *) fx->extract("a", dim);
    double *pb = (double *) fx->extract("b", dim);
    double *pf = (double *) fx->extract("f", dim);
    if (!pa || !pb || !pf) continue;

    double zl, zt, zh;
    if (!landscape_extrema(*pa, *pb, *pf, zl, zt, zh)) continue;

    // Theta at the two minima, evaluated at the crossover this landscape actually gets
    double zc = zstar_auto ? zt : z_star;
    double th_lo = 1.0 / (1.0 + exp(-sigmoid_alpha * (zl - zc)));
    double th_hi = 1.0 / (1.0 + exp(-sigmoid_alpha * (zh - zc)));
    double res = lam_max_all * th_lo * th_lo;
    double phi_lo = *pa * zl * zl * zl * zl - *pb * zl * zl + *pf * zl;
    double phi_top = *pa * zt * zt * zt * zt - *pb * zt * zt + *pf * zt;
    double disp = fabs(zt - zc);

    worst_res = MAX(worst_res, res);
    worst_disp = MAX(worst_disp, disp * sigmoid_alpha);

    if (comm->me == 0)
      utils::logmesg(lmp,
                     "  fix {}: a={} b={} f={}\n"
                     "    inactive min z- = {:.6}  Theta = {:.4e}\n"
                     "    barrier top  z0 = {:.6}  crossover at {:.6}  |z0-zc| = {:.4e}\n"
                     "    bonding  min z+ = {:.6}  1-Theta = {:.4e}\n"
                     "    barrier f+ = {:.6} kT   residual inactive-inactive lambda = {:.4e}\n",
                     fx->id, *pa, *pb, *pf, zl, th_lo, zt, zc, disp, zh, 1.0 - th_hi,
                     phi_top - phi_lo, res);
  }

  if (worst_res > 100.0 * theta_tol)
    error->all(FLERR,
               "Pair style lj/cut/mlo: the residual attraction between two inactive "
               "particles is {} (tolerance {}), so inactive states are not effectively "
               "non-attracting. Increase alpha, or raise theta_check",
               worst_res, theta_tol);
  if (worst_res > theta_tol && comm->me == 0)
    error->warning(FLERR,
                   "Pair style lj/cut/mlo: the residual attraction between two inactive "
                   "particles is {}, above the tolerance {}. Increase alpha to suppress it",
                   worst_res, theta_tol);
  if (worst_disp > 1.0 && comm->me == 0 && !zstar_auto)
    error->warning(FLERR,
                   "Pair style lj/cut/mlo: the Theta crossover is displaced from the barrier "
                   "top by {} times its own width 1/alpha. The interaction's notion of the "
                   "chemical state disagrees with the landscape's. Use `zstar auto`",
                   worst_disp);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double MLOPairLJCut::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i], epsilon[j][j], sigma[i][i], sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i], sigma[j][j]);
    // The in-plane cutoff mixes like a distance, matching how stock lj/cut mixes `cut`.
    cut_2d[i][j] = mix_distance(cut_2d[i][i], cut_2d[j][j]);
    // lambda_max mixes GEOMETRICALLY by default, not with the arithmetic rule the
    // published Ashbaugh-Hatch/HPS models use. Two reasons. It is self-consistent with the
    // product form lambda_ij = Theta_i*Theta_j already chosen here, so lambda_max becomes a
    // per-type "bonding valence": lambda_ij = (sqrt(L_i) Theta_i)(sqrt(L_j) Theta_j). And
    // it makes lambda_ii = 0 mean inert with EVERY partner under wildcard mixing, whereas
    // the arithmetic rule would give lambda_12 = 0.5 for an inert type paired with a fully
    // attractive one -- silently reintroducing attraction to a species declared inert.
    if (lambda_mix_flag == 0)
      lambda_max[i][j] = sqrt(lambda_max[i][i] * lambda_max[j][j]);
    else
      lambda_max[i][j] = 0.5 * (lambda_max[i][i] + lambda_max[j][j]);
  }

  if (cut_2d[i][j] <= 0.0)
    error->all(FLERR, "Pair style lj/cut/mlo has no in-plane cutoff for types {} {}", i, j);

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);

  double CUBE_ROOT_TWO = pow(2.0, 1.0 / 3.0);
  LJ_MINIMUM[i][j] =
      CUBE_ROOT_TWO * sigma[i][j] * sigma[i][j];    // Square Minimum of the LJ potential in 2D

  cut_2d_sq[i][j] = cut_2d[i][j] * cut_2d[i][j];

  // A cutoff inside the LJ minimum truncates the potential on the repulsive wall, so the
  // force jumps by a large amount at r_c. Almost certainly a mistake; warn rather than
  // forbid, since a deliberately tiny purely-repulsive pair is conceivable.
  if (cut_2d_sq[i][j] < LJ_MINIMUM[i][j] && comm->me == 0)
    error->warning(FLERR,
                   "Pair style lj/cut/mlo in-plane cutoff {} for types {} {} is inside the "
                   "LJ minimum {}, so the potential is truncated on its repulsive wall",
                   cut_2d[i][j], i, j, sqrt(LJ_MINIMUM[i][j]));

  // FIXED by CLAUDE (2026-08-31): the energy shift has to refer to the radius at which the
  // potential is actually truncated, which is the in-plane cutoff - not cut[i][j], which
  // only sets the (deliberately inflated) neighbor cutoff. Now per type pair.
  if (offset_flag) {
    double ratio = sigma[i][j] / cut_2d[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio, 12.0) - pow(ratio, 6.0));
  } else
    offset[i][j] = 0.0;

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];
  // FIXED by CLAUDE (2026-08-31): compute() reads epsilon and LJ_MINIMUM as
  // [itype][jtype] with no ordering guarantee, but only the upper triangle was ever
  // filled - so for more than one atom type the lower triangle was read uninitialized.
  epsilon[j][i] = epsilon[i][j];
  sigma[j][i] = sigma[i][j];
  LJ_MINIMUM[j][i] = LJ_MINIMUM[i][j];
  cut_2d[j][i] = cut_2d[i][j];
  cut_2d_sq[j][i] = cut_2d_sq[i][j];
  // Mirrored because compute() indexes [itype][jtype] with no ordering guarantee, and
  // because fix adapt writes only the upper triangle before calling Pair::reinit(), which
  // re-runs init_one() for j >= i only.
  lambda_max[j][i] = lambda_max[i][j];

  // Long-range tail corrections are rejected in init_style(): the 3d formula does not
  // apply to an interaction that acts in the xy plane only and is scaled by lambda(z).

  // The returned value is the NEIGHBOR cutoff, which Pair::init() turns into both
  // cutsq[i][j] and cutforce. It has to be inflated over the in-plane cutoff by the full
  // z budget, otherwise an in-plane neighbour separated in z drops out of the list.
  //
  // Bit-exactness shim: when the budget was derived from the legacy positional arguments
  // and this pair uses the default in-plane cutoff, return cut_global verbatim rather than
  // its sqrt round-trip. The two agree to ~1 ulp, but a last-bit difference could flip a
  // neighbour sitting exactly at the cutoff, and existing results should reproduce exactly.
  if (!zbudget_explicit && cut_2d[i][j] == cutoff_2d)
    cut[i][j] = cut_global;
  else
    cut[i][j] = sqrt(cut_2d[i][j] * cut_2d[i][j] + z_budget * z_budget);

  cut[j][i] = cut[i][j];

  // Internal assertion: every pair inside the in-plane cutoff must survive the 3d neighbor
  // build. True by construction above, but cheap to guard against a future edit to the
  // derivation. The tolerance absorbs the shim's 1-ulp disagreement with the sqrt.
  if (cut[i][j] * (1.0 + 1.0e-12) < sqrt(cut_2d[i][j] * cut_2d[i][j] + z_budget * z_budget))
    error->all(FLERR,
               "Pair style lj/cut/mlo internal error: derived neighbor cutoff {} for types "
               "{} {} is too small for an xy cutoff of {} and a z budget of {}",
               cut[i][j], i, j, cut_2d[i][j], z_budget);

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
        // cut_2d, the in-plane interaction range -- NOT cut[i][j], which is derived from
        // it and the z budget and would be meaningless to restore directly.
        fwrite(&cut_2d[i][j], sizeof(double), 1, fp);
        fwrite(&lambda_max[i][j], sizeof(double), 1, fp);
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
          if (restart_version >= 2) {
            utils::sfread(FLERR, &cut_2d[i][j], sizeof(double), 1, fp, nullptr, error);
            utils::sfread(FLERR, &lambda_max[i][j], sizeof(double), 1, fp, nullptr, error);
          } else {
            // Pre-versioning files stored cut[i][j], the deliberately inflated NEIGHBOR
            // cutoff. Restoring it as an in-plane cutoff would silently reproduce the
            // 5.59-sigma interaction range of CODE_REVIEW P1 #1, so read and discard it.
            double legacy_cut;
            utils::sfread(FLERR, &legacy_cut, sizeof(double), 1, fp, nullptr, error);
            cut_2d[i][j] = -1.0;
            lambda_max[i][j] = 1.0;
          }
        }
        MPI_Bcast(&epsilon[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&sigma[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut_2d[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&lambda_max[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void MLOPairLJCut::write_restart_settings(FILE *fp)
{
  // All pair_style settings are stored now, so read_restart + run 0 reproduces the run
  // without a re-issued pair_style. A leading 64-bit magic makes pre-versioning files
  // detectable on read: their first 8 bytes are the double cut_global, and a fixed 64-bit
  // pattern collides with that only with probability 2^-64. (A 32-bit magic would risk
  // colliding with the low mantissa word.)
  const uint64_t magic = MLO_RESTART_MAGIC;
  const int version = 2;
  fwrite(&magic, sizeof(uint64_t), 1, fp);
  fwrite(&version, sizeof(int), 1, fp);

  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&z_star, sizeof(double), 1, fp);
  fwrite(&cutoff_2d, sizeof(double), 1, fp);
  fwrite(&sigmoid_alpha, sizeof(double), 1, fp);
  fwrite(&z_budget, sizeof(double), 1, fp);
  fwrite(&theta_tol, sizeof(double), 1, fp);
  fwrite(&zbudget_explicit, sizeof(int), 1, fp);
  fwrite(&lambda_mix_flag, sizeof(int), 1, fp);
  fwrite(&zstar_auto, sizeof(int), 1, fp);
  fwrite(&zcheck_flag, sizeof(int), 1, fp);
  fwrite(&zperiodic_ok, sizeof(int), 1, fp);
  fwrite(&theta_check_on, sizeof(int), 1, fp);
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
  restart_version = 1;

  if (me == 0) {
    long pos = ftell(fp);
    uint64_t magic = 0;
    utils::sfread(FLERR, &magic, sizeof(uint64_t), 1, fp, nullptr, error);
    if (magic == MLO_RESTART_MAGIC) {
      utils::sfread(FLERR, &restart_version, sizeof(int), 1, fp, nullptr, error);
      if (restart_version != 2)
        error->one(FLERR,
                   "Pair style lj/cut/mlo restart file has unsupported format version {}",
                   restart_version);
      utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
      utils::sfread(FLERR, &z_star, sizeof(double), 1, fp, nullptr, error);
      utils::sfread(FLERR, &cutoff_2d, sizeof(double), 1, fp, nullptr, error);
      utils::sfread(FLERR, &sigmoid_alpha, sizeof(double), 1, fp, nullptr, error);
      utils::sfread(FLERR, &z_budget, sizeof(double), 1, fp, nullptr, error);
      utils::sfread(FLERR, &theta_tol, sizeof(double), 1, fp, nullptr, error);
      utils::sfread(FLERR, &zbudget_explicit, sizeof(int), 1, fp, nullptr, error);
      utils::sfread(FLERR, &lambda_mix_flag, sizeof(int), 1, fp, nullptr, error);
      utils::sfread(FLERR, &zstar_auto, sizeof(int), 1, fp, nullptr, error);
      utils::sfread(FLERR, &zcheck_flag, sizeof(int), 1, fp, nullptr, error);
      utils::sfread(FLERR, &zperiodic_ok, sizeof(int), 1, fp, nullptr, error);
      utils::sfread(FLERR, &theta_check_on, sizeof(int), 1, fp, nullptr, error);
    } else {
      // pre-versioning layout: rewind and read it as it was written
      restart_version = 1;
      fseek(fp, pos, SEEK_SET);
      utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    }
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tail_flag, sizeof(int), 1, fp, nullptr, error);
  }

  MPI_Bcast(&restart_version, 1, MPI_INT, 0, world);
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&tail_flag, 1, MPI_INT, 0, world);

  if (restart_version >= 2) {
    MPI_Bcast(&z_star, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&cutoff_2d, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&sigmoid_alpha, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&z_budget, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&theta_tol, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&zbudget_explicit, 1, MPI_INT, 0, world);
    MPI_Bcast(&lambda_mix_flag, 1, MPI_INT, 0, world);
    MPI_Bcast(&zstar_auto, 1, MPI_INT, 0, world);
    MPI_Bcast(&zcheck_flag, 1, MPI_INT, 0, world);
    MPI_Bcast(&zperiodic_ok, 1, MPI_INT, 0, world);
    MPI_Bcast(&theta_check_on, 1, MPI_INT, 0, world);
    // A v2 file carries everything, so no re-issued pair_style is required.
    settings_set = 1;
  } else if (comm->me == 0) {
    error->warning(FLERR,
                   "Pair style lj/cut/mlo is reading a pre-versioning restart file, which "
                   "does not store z_star, the xy cutoff, alpha or the z budget. Reissue "
                   "the pair_style and pair_coeff commands");
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void MLOPairLJCut::write_data(FILE *fp)
{
  // Keyword form, so read_data can feed these straight back to coeff(), which rejects a
  // bare 5th positional precisely because the old format put the neighbour cutoff there.
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp, "%d %g %g cut %g lambda %g\n", i, epsilon[i][i], sigma[i][i], cut_2d[i][i],
            lambda_max[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void MLOPairLJCut::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g cut %g lambda %g\n", i, j, epsilon[i][j], sigma[i][j],
              cut_2d[i][j], lambda_max[i][j]);
}

/* ---------------------------------------------------------------------- */

double MLOPairLJCut::single(int /*i*/, int /*j*/, int /*itype*/, int /*jtype*/, double /*rsq*/,
                            double /*factor_coul*/, double /*factor_lj*/, double &/*fforce*/)
{
  // FIXED by CLAUDE (2026-08-31): this used to return a hardcoded 10.0 while
  // single_enable was left at its default of 1, so compute group/group, fix bond/create
  // or a fix widom/gcmc line without `full_energy` silently got a plausible wrong number.
  // single() cannot be implemented from rsq alone: the energy depends on the in-plane
  // separation and on z_i and z_j separately, not on the 3d pair distance.
  error->all(FLERR,
             "Pair style lj/cut/mlo does not support single(). Its energy is not a "
             "function of the 3d pair distance. Use the full_energy option where available");
  return 0.0;
}

/* ---------------------------------------------------------------------- */

void MLOPairLJCut::born_matrix(int /*i*/, int /*j*/, int /*itype*/, int /*jtype*/,
                               double /*rsq*/, double /*factor_coul*/, double /*factor_lj*/,
                               double &/*dupair*/, double &/*du2pair*/)
{
  // FIXED by CLAUDE (2026-08-31): the inherited body evaluated the plain 3d LJ
  // derivatives, ignoring both the xy-only distance and the lambda(z) prefactor.
  // born_matrix_enable is now 0; error out rather than return wrong numbers.
  error->all(FLERR, "Pair style lj/cut/mlo does not support born_matrix()");
}

/* ---------------------------------------------------------------------- */

void *MLOPairLJCut::extract(const char *str, int &dim)
{
  // dim has to be set per key, not once up front: the matrices are dim 2 and the globals
  // are dim 0.
  dim = 2;
  if (strcmp(str, "epsilon") == 0) return (void *) epsilon;
  if (strcmp(str, "sigma") == 0) return (void *) sigma;
  // The attraction matrix is the natural knob for a non-equilibrium ramp via fix adapt.
  // Two caveats, both inherent to the fix adapt interface: it writes only the upper
  // triangle (init_one mirrors it, so that is handled), and Pair::reinit() re-mixes any
  // pair with setflag == 0, so adapting a MIXED pair is silently a no-op - adapt the
  // diagonals and let the mixing rule propagate, or set every pair with pair_coeff.
  if (strcmp(str, "lambda_max") == 0) return (void *) lambda_max;
  if (strcmp(str, "lambda") == 0) return (void *) lambda_max;

  dim = 0;
  if (strcmp(str, "z_star") == 0) return (void *) &z_star;
  if (strcmp(str, "alpha") == 0) return (void *) &sigmoid_alpha;
  if (strcmp(str, "z_budget") == 0) return (void *) &z_budget;

  // cut_2d and cut are deliberately NOT exposed: Pair::reinit() calls init_one() but never
  // re-runs Neighbor::init(), so adapting an interaction range would leave cutneighsq
  // stale and silently drop pairs.
  return nullptr;
}

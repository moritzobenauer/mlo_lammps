/* -*- c++ -*- ----------------------------------------------------------
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
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/mlo,MLOPairLJCut)
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_MLO_H
#define LMP_PAIR_LJ_CUT_MLO_H

#include "pair.h"

namespace LAMMPS_NS {

class MLOPairLJCut : public Pair {
 public:
  MLOPairLJCut(class LAMMPS *);
  ~MLOPairLJCut() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void born_matrix(int, int, int, int, double, double, double, double &, double &) override;
  void *extract(const char *, int &) override;

//   void compute_inner() override;
//   void compute_middle() override;
//   void compute_outer(int, int) override;

 protected:
  double cut_global;    // arg[0]: legacy 3D neighbor cutoff; only seeds z_budget now
  double z_star;    // z position of the free energy barrier, where Theta(z) transitions from 0 to 1
  double cutoff_2d; // arg[2]: DEFAULT in-plane interaction range; seeds cut_2d[i][j]
  double sigmoid_alpha;    // steepness of Theta(z); optional 4th pair_style argument
  // Guaranteed |dz| headroom. The interaction is gated on the in-plane distance while the
  // neighbor list is built from true 3D distances, so a pair close in xy but far apart in z
  // drops out of the list and loses even its repulsive core. Every pair's neighbor cutoff is
  // therefore sqrt(cut_2d^2 + z_budget^2), and the z spread is checked against z_budget.
  double z_budget;
  int zbudget_explicit;    // 1 if the zbudget keyword was given
  int settings_set;        // 1 once settings() has run (replaces the cutoff_2d < 0 sentinel)
  int lambda_mix_flag;     // 0 = geometric (default), 1 = arithmetic
  int zstar_auto;          // 1 = derive the Theta crossover from fix free_energy (default)
  int zcheck_flag;         // 0 = off, 1 = warn (default), 2 = strict
  int zperiodic_ok;        // user acknowledgement of the periodic-z hazards
  int theta_check_on;      // 1 = run the inactive-residual diagnostic
  double theta_tol;        // tolerance on the residual inactive-inactive lambda
  bigint last_zcheck;      // dedupe the per-rebuild z guard
  int zwarn_issued;        // throttle the rank-local warning to once per run
  int restart_version;     // 1 = pre-versioning restart file, 2 = current
  double **cut;            // DERIVED neighbor cutoff = sqrt(cut_2d^2 + z_budget^2)
  double **cut_2d, **cut_2d_sq;    // per-type-pair in-plane interaction range
  // Per-type-pair attraction scale. epsilon owns the repulsive core, lambda_max owns the
  // attraction, so lambda_max = 0 gives exact WCA at full epsilon.
  double **lambda_max;
  // Per-ATOM-TYPE crossover of Theta(z). The landscape phi(z) is already per type in
  // practice (two_particles.in runs two fix free_energy instances on type-based groups),
  // so a single global z_star cannot sit at the barrier top of both.
  double *z_star_type;
  double **epsilon, **sigma;
  double **lj1, **lj2, **lj3, **lj4, **offset, **LJ_MINIMUM;

  virtual void allocate();

  // Guards and diagnostics (see the implementations for the reasoning).
  void check_z_budget();
  void derive_z_star();
  void report_theta_residual();
  // Solves phi'(z) = 4a z^3 - 2b z + f = 0 by bisection. Returns 1 and fills the three
  // extrema when phi is a genuine double well, 0 otherwise.
  int landscape_extrema(double a, double b, double f, double &zlo, double &ztop, double &zhi);
};

}    // namespace LAMMPS_NS

#endif
#endif

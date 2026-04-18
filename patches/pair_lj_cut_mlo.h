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
  double cut_global;
  double z_star;    // z position of the free energy barrier, where Theta(z) transitions from 0 to 1
  double cutoff_2d; // cutoff for the 2D distance in the xy plane
  double **cut;
  double **epsilon, **sigma;
  double **lj1, **lj2, **lj3, **lj4, **offset, **LJ_MINIMUM;
  double *cut_respa;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

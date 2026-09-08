#ifdef FIX_CLASS
// This line MUST match the name you use in the input script
FixStyle(free_energy, FixFreeEnergy); 
#else

#ifndef LMP_FIX_FREEENERGY_H
#define LMP_FIX_FREEENERGY_H

#include "fix.h"

namespace LAMMPS_NS {


class FixFreeEnergy : public Fix {
 public:
  FixFreeEnergy(class LAMMPS *, int, char **);
  ~FixFreeEnergy() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_integrate() override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  // FIXED by CLAUDE (2026-08-31): phi(z) is now also applied during minimization,
  // and the frozen-z DOF are reported to temperature computes.
  void min_post_force(int) override;
  bigint dof(int) override;
  double compute_scalar() override;

  // Exposes the landscape coefficients of phi(z) = a*z^4 - b*z^2 + f*z so that
  // pair_style lj/cut/mlo can locate the barrier top and tie the Theta(z) crossover
  // to the actual free-energy landscape instead of a hand-supplied z_star.
  void *extract(const char *, int &) override;

 protected:
  void freeze_z_velocities();

  double coeff_a, coeff_b, coeff_f; // coefficients for z^4, z^2, z
  int disable_reactions;
  int ilevel_respa;
  int eflag;
  double e_total, e_total_all;
};

}

#endif
#endif

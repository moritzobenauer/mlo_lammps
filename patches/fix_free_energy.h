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
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  double compute_scalar() override;

 protected:
  double coeff_a, coeff_b, coeff_f; // coefficients for z^4, z^2, z
  int disable_reactions;
  int ilevel_respa;
  int eflag;
  double e_total, e_total_all;
};

}

#endif
#endif

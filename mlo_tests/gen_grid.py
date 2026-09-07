#!/usr/bin/env python3
"""Regenerate grid.data: a 14x14 xy grid with random z.

The grid matters. `create_atoms random ... overlap` rejects overlaps in 3D, which does NOT
guarantee in-plane separation for this pair style -- two atoms can be 0.6 apart in 3D while
nearly coincident in xy, which is a divergent in-plane force. Placing atoms on an xy grid
makes the minimum in-plane separation exact.
"""
import random
random.seed(20260907)
n, a = 14, 1.35
L = n * a
atoms = []
for i in range(n):
    for j in range(n):
        t = 1 if (i * n + j) % 5 else 2
        atoms.append((t, i*a + 0.15*(random.random()-0.5),
                         j*a + 0.15*(random.random()-0.5),
                         random.uniform(-1.0, 1.0)))
with open("grid.data", "w") as f:
    f.write("mlo mpi consistency test\n\n%d atoms\n2 atom types\n\n" % len(atoms))
    f.write("0.0 %.10f xlo xhi\n0.0 %.10f ylo yhi\n-6.0 6.0 zlo zhi\n\n" % (L, L))
    f.write("Masses\n\n1 1.0\n2 1.0\n\nAtoms # atomic\n\n")
    for k, (t, x, y, z) in enumerate(atoms, 1):
        f.write("%d %d %.12f %.12f %.12f\n" % (k, t, x, y, z))

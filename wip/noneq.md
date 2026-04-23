# fix that adds a constant force along z to mimic an external reservoir

To-Do and Idea
- () code `fix external_chemostat()` which takes an atom group $G$ and applies a constant force along $z$: `f[i][2] += Delta` where $Delta$ is a measure for the non-equilibrium extent of the system. For $\Delta  = 0.0$ the system should reach thermodynamic equilibrium!
- () external chemostat should be a continously decreasing function of $z$: maybe $f(z) = \Delta \exp(-z)$
- () pair-wise interactions need to be able to account for pair-wise chemical reactions in the following form
  $$ f_z(r)^\mathrm{monomer} = \frac{\zeta_{ij} \lambda_{ij}}{1+\exp(\alpha(r-r_c))} $$

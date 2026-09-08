1. There is a `patches/patches` directory which is wrong. Everything should be in the parent directory, otherwise the della patching is not working.
2. Making sure that input scripts do not have periodic z boundary conditions. For a finite phi(z) potential we should never leave the physically reasonable domain anyway. But better safe than sorry.
3. The claude modified cpp files seem a bit cluttered. I need to clean up.
4. Local catalysis implementation
5. Chemical driving with detailed balance???


Analysis:

[] Automate the jump analysis after a succesful run --> report JUMPS, jump rate, active to inactive, prob_active, prob_inactive, ratio_active_to_inactive, ratio_kinetic_active_to_inactive --> This is a good indicator to see if we have reached steady state as well!
[] In the long term, the distribution P(j) is probably interesting. I do not know how, yet. Is it P(j(r))? To correlate local density with jump probability. But this should probably be done by a local posterior script.
[] automatically extract P(z)! That is my PMF and should resemble phi(z) in the absence of any attractive interactions! If there are attractive interactions, then P(z) is a function of the interaction parameter and density. But I need to quanitfy that.


The big question:

Quantifying the difference between <J(P(z))> and J(<P(z)>). When is it a good approximation to get P(z) and then let a 1D particle hop along this potential?



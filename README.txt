TL FITTING
===================

This folder contains the functions tlfit.m and rlfit.m to calculate the fitting parameters 
SL (source level), GLF (geometric loss factor) and ALF (attenuation loss factor) of the
curves TL = GLF*log10(R) + ALF*R and RL = SL - GLF*log10(R) - ALF*R, using data from 
measured or simulated range-dependent transmission loss TL(R) and received level RL(R) curves.

The functions rely on a box-constrained version of fminsearch.m, named fminsearchbnd.m. The 
latter is inspired in a function of the same name, but be aware that tlfit.m and rlfit.m will 
only work with my version.

Tests for tlfit.m and rlfit.m are included in folder "Tests (.m)".

Other functions with potential future use are included in folder "Other (.m)".

Documents that were relevant for the creation of fminsearchbnd.m, addressing the conversion
between unconstrained and constrained optimisation problems, are included in folder "Docs".

Old versions of current functions are included in folder "Archive (.m)".
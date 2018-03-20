# FHDq2D
Fluctuating hydrodynamics of binary mixtures in quasi2D geometries

This MATLAB code solves the equations of fluctuating hydrodynamics for a binary mixture of dynamically-identical (labeled or tagged) spherical colloids in quasi-2D confinement. It can solve the nonlinear FHD Eq. (53) in the paper:

[**Hydrodynamic fluctuations in quasi-two dimensional diffusion**](https://cims.nyu.edu/~donev/FluctHydro/Quasi2D_DDFT.pdf), R. P. Peláez, F. Balboa Usabiaga, S. Panzuela, Q. Xiao, R. Delgado-Buscalioni and A. Donev, [ArXiv:1802.07356](https://arxiv.org/abs/1802.07356), submitted to J. Stat. Mech., 2018.

but it can also solve a linearized version of this equation.

The hydrodynamic kernel **R** can be changed via the flag velMode (1 for quasi2D, 2 for true2D velocity type, or 3 for Saffman membrane hydrodynamics) in main.m. As explained in the paper, for quasi2D we really solve Eq. (52) in the paper, instead of (53).

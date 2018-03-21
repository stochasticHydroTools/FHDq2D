# FHDq2D
Fluctuating hydrodynamics of binary mixtures in quasi2D geometries

This MATLAB code solves the equations of fluctuating hydrodynamics for a binary mixture of dynamically-identical (labeled or tagged) spherical colloids in quasi-2D confinement. It can solve the nonlinear FHD equations proposed in the paper:

[**Hydrodynamic fluctuations in quasi-two dimensional diffusion**](https://cims.nyu.edu/~donev/FluctHydro/Quasi2D_DDFT.pdf), R. P. Peláez, F. Balboa Usabiaga, S. Panzuela, Q. Xiao, R. Delgado-Buscalioni and A. Donev, [ArXiv:1802.07356](https://arxiv.org/abs/1802.07356), submitted to J. Stat. Mech., 2018.

The equations are summarized in the documentation file GiantFluct_2D.pdf. The full nonlinear equation is Eq. (3) in this documentation. But the code can also solve a linearized version of this equation, see Eq. (13) in the documentation. As explained in the paper, for quasi2D we really solve the simplified Eq. (4) in the documentation, not Eq. (3).

The hydrodynamic kernel **R** can be changed via the flag velMode in main.m (1 for quasi2D, -1 for quasi2D with incompressible component only, 2 for true2D velocity type, or 3 for (2+1)d Saffman hydrodynamics). 

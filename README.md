![hyper2D-banner](https://raw.githubusercontent.com/wiki/sbocce/hyper2d/imgs/hyper2D_banner.png)

Hyper2D is a minimalistic finite volume solver for hyperbolic PDEs and hypersonic flows.

Check out the Wiki: 
https://github.com/sbocce/hyper2d/wiki

### Citing Hyper2D in your work

Hyper2D comes for free. 
We don't have a specific journal publication for it, however, if you use Hyper2D in your research, 
you could cite the following publication, where we use Hyper2D to solve the Euler equations and some other
hyperbolic PDEs:

S. Boccelli, J. G. McDonald, T. E. Magin, _14-moment maximum-entropy modelling of collisionless ions for Hall thruster discharges_, http://arxiv.org/abs/2202.04159

### Quick intro

Hyper2D is designed to be as simple as possible and suitable for teaching purposes.
You will find in this package various versions of Hyper2D, with increasing capabilities and complexity.
Each version is completely independent, and is compiled and run independently from the others.

### Requirements

Hyper2D was tested on GNU/Linux. The software requirements depend on the specific Hyper2D 
implementation to be used:

1. hyper1D directory:

  - PDE1D_Octave -> runs under Octave and/or MATLAB (with little or no modifications)
  - src_hyper1D and src_hyper1D_simple -> Fortran versions, require gfortran and make.

2. hyper2D_single_core directory: Fortran versions, require gfortran and make

3. hyper2D_CUDA: requires a GPU. Software requirements: make and some CUDA implementation 
  (tested working with the NVidia HPC SDK v22, that include the nvfortran compiler;  
  other CUDA compilers may work, but may require to modify the parameters in the Makefile).

### Compiling and running

First, you can clone this repository from github with:

$ git clone github.com/sbocce/hyper2D

This package includes various versions of Hyper2D, with increasing complexity. 
Once you have identified the version that you are interested in (1D, 2D single core or 
2D CUDA), go to the desired src directory: you will find a Makefile. Compile with:

$ make

This will create an executable file "hyper2D.exe" in the local directory.
Hyper2D will print its output in a local "dumps" directory. Make sure that it
exists: 

$ mkdir -p dumps

You are now ready to run the program:

$ ./hyper2D.exe

After applying modifications to the program, such as changing the simulated conditions,
you need to compile again.
Suggestion:
the following command allows to clean the "dumps" directory, removing the output
of previous simulations; then compile the sources and run hyper2D:

$ make cleanoutput && make && ./hyper2D.exe

### Capabilities and versions of Hyper2D

As mentioned, Hyper2D comes in different implementations, with increasing degrees of complexity.
Here is a summary. 

The **hyper1D directory** contains the following:

* PDE1D

  This is a simple one-dimensional finite volume solver, implemented in Octave/MATLAB,
  aimed at familiarizing with the Finite Volume Method. It solves the Euler equations 
  for the Sod shock tube problem.

* hyper1D_simple

  The simplest Fortran implementation. Pretty much a direct re-implementation of the 
  Octave/MATLAB version.

* hyper1D

  Fortran version of PDE1D, including higher order (in space and time) numerics.

------ The hyper2D_single_core contains simple 2D implementations that run on a single core ------

#### hyper2D_basic

  This is the simplest implementation of Hyper2D, in Fortran. It includes the Shallow 
  Water equations and the Euler equations (either with constant specific heat or with 
  molecular vibration for diatomic molecules). You will find these files in the 
  "various_PDEs" directory. You just have to rename the files  pde_shallow_water.f03 
  or pde_euler.f03 to "pde.f03", put it in the main "src" directory and compile.
 
  This version of Hyper2D solves the problem in a Cartesian grid, using a first order 
  method  in space and explicit in time. A spatially second order version, using 
  van Leer's MUSCL approach and TVD slope limiter is available in a sub-directory.
 
  The domain is empty, and one can study the effect of boundaries or initialize the 
  solution with a non-trivial initial state and track its evolution.

  Includes an axisymmetric formulation.

#### hyper2D_flatplate

  This version builds on hyper2D_basic and introduces a flat plate inside the domain.
  It uses a Cartesian grid, and you can easily introduce other flat surfaces that
  are aligned with grid cells. 
 
#### hyper2D_hypersonic

  This directory contains a non-equilibrium two-temperatures solver for the Euler 
  equations, targeted at hypersonic simulations. Two temperatures are considered: 
  a roto-translational temperature T and a vibrational temperature Tv. 
  These temperatures relax following the Landau-Teller model.
  
  The time scales for vibrational relaxation can be very fast. Therefore, implicit time
  integration could be preferred in some conditions. The directory "other_files" contains
  both an explicit and a point-implicit integration scheme. The point-implicit scheme
  is only done on the vibrational energy equation, since the source term for the other
  equations is zero.

#### hyper2D_kinetic
  
  This version solves the kinetic Boltzmann/Vlasov equation, in a "1D in physical space" 
  and "1V in velocity" phase space. The variable "x" represents physical space and 
  the variable "y" represents velocity.

#### hyper2D_unstructured
  
  This version reads an unstructured grid of quadrilaterals in SU2 format. The grid
  is pre-processed using some Octave/MATLAB tools available in the "mesh_preprocessing"
  directory, resulting in a "mesh.hyp" and a "nodes.hyp" files that are intelligible 
  by Hyper2D.

----- The hyper2D_CUDA directory: -----------

#### hyper2D_CUDA_simplest

  The simplest GPU implementation of Hyper2D, written in CUDA Fortran, intended for 
  familiarizing with CUDA instructions. Runs on NVidia GPUs, on a single GPU.

#### hyper2D_CUDA
 
  This is the most advanced version of hyper2D, is written in CUDA Fortran and allows 
  you to run Hyper2D on an NVidia GPU, with higher order numerics. 
  The code runs on a single GPU. The implementation includes axial symmetry (see the 
  README in the specific folder).



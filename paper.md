---
title: 'Hyper2D: a finite-volume solver for hyperbolic PDEs and hypersonic flows'
tags:
  - hyperbolic PDEs
  - hypersonic flows
  - finite-volume method
authors:
  - name: Stefano Boccelli^[corresponding author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0003-2276-9606
    affiliation: "1, 2, 3" # (Multiple affiliations must be quoted)
  - name: James G. McDonald # note this makes a footnote saying 'co-first author'
    affiliation: 1
  - name: Thierry E. Magin
    affiliation: 3
affiliations:
 - name: University of Ottawa, ON, Canada # Lyman Spitzer, Jr. Fellow, Princeton University
   index: 1
 - name: Politecnico di Milano, Italy
   index: 2
 - name: von Karman Institute for Fluid Dynamics, Belgium
   index: 3
date: 6 January 2022
bibliography: paper.bib

# NO ## # Optional fields if submitting to a AAS journal too, see this blog post:
# NO ## # https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# NO ## aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# NO ## aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Conservation equations from various branches of mathematical physics 
often take the form of hyperbolic partial differential equations (PDEs).
The finite-volume method (FVM), see @leveque2002finite, allows to obtain an approximate solution, by 
first discretizing these equations on a spatial grid, and then marching in time in order to 
follow the time evolution of the discretized fields.

Many different conservation equations can be written as a set of hyperbolic PDEs.
Some examples are the scalar Burgers equation, the shallow water equations,
the Euler equations for gas dynamics and countless other systems.

# Statement of need

`Hyper2D` is a minimalistic FVM solver for hyperbolic equations, written in Fortran 2003.
Its simple structure makes `Hyper2D` accessible and easy to learn.
In particular, `Hyper2D` targets 

1. students and newcomers to the finite-volume method;

2. researchers wishing to implement their own set of PDEs.

### An introductory code

The `Hyper2D` package contains a number of independent implementations of the FVM solver.
Each standalone implementation has an increasing level of complexity and offers further features.
Students can step in at the appropriate level, according to their knowledge, and work their way
towards the more thorough implementations.

The simplest version is one-dimensional, written in Octave/MATLAB, and allows the students
to familiarize with the FVM method.
Two-dimensional versions are then provided, for simple Cartesian grids, as well as unstructured
grids of quadrilateral elements.
In each version, additional complexities are introduced, such as internal solid surfaces, more 
advanced physics or numerics.
Finally, a CUDA Fortran version is provided, that allows to exploit GPGPU acceleration.

Currently, `Hyper2D` is being used during the Fundamentals of Hypersonic Flows lessons at 
Politecnico di Milano.
Students typically use it for testing out the effect of the thermochemical model on the flow topology, 
for implementing new physical models and learning the numerics.

### A research code

The simple structure of `Hyper2D` makes it easily customizable.
One can easily implement a new set of PDEs, by modifying the simple examples provided.
Object oriented programming is avoided, and the code is completely procedural.
The whole solution vector is available throughout most of the code:
this gives additional flexibility, as one has the freedom to implement non-local source 
terms for instance.

Until now, `Hyper2D` has been employed in the following publications, for solving moment equations for non-equilibrium 
fluid dynamics, @boccelli2021thesis, and plasma physics problems related to space propulsion 
@boccelli2022maxention, and @boccelli2020collisionless.

# Governing equations and numerics

`Hyper2D` solves governing equations in the form

\begin{equation}\label{eq:PDE}
  \frac{\partial U}{\partial t}
+ \frac{\partial F_x}{\partial x}
+ \frac{\partial F_y}{\partial y}
=
  S \, ,
\end{equation}

where $U$ is the state, $F_x$ and $F_y$ are the fluxes in the $x$ and $y$ directions and $S$ are source terms. 
These quantities are scalars (in case of a scalar equatios) or vectors (for a system of PDEs).

`Hyper2D` offers some different numerical schemes for discretizing and solving \autoref{eq:PDE}. 
Such schemes are independent on the specific PDEs to be solved, and should work out of the box even after 
new equations are implemented.
Currently, `Hyper2D` offers the following explicit time integrators:

- 1st order Forward Euler;

- 2nd order Midpoint Euler;

- 3rd order Runge-Kutta;

- 3rd order low-storage TVD Runge-Kutta.

And the computation of convective fluxes has the following options: 

- 1st order in space, no reconstruction;

- 2nd order van Leer's MUSCL with TVD slope limiting;

- 3rd order WENO;

- 5th order WENO.


# Example

As a practical example, \autoref{fig:jet} shows the simulation of a supersonic axisymmetric jet,
obtained solving the relativistic Euler equations.
The initial conditions are taken from the test case by @mignone2005hllc.

![Simulation of a relativistic axisymmetric jet. Logarithm of the density in the rest frame.\label{fig:jet}](jet1.png)

# Acknowledgements

The authors would like to thank Prof. Aldo Frezzotti (Politecnico di Milano) for 
the insightful discussions about the hypersonic non-equilibrium features implemented 
in Hyper2D.

# References
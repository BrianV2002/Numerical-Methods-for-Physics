# Burgers Equation in Fortran (1D)

This repository contains a **numerical implementation** of the 1D Burgers
equation using Fortran The goal is that any reader can
understand what is being solved, why it works, and how to run the
code without having to dig into every line of Fortran.

The main files are:

- `burgers.f90`  → main program (time integration of the equation)
- `mod.f90`      → module with parameters and auxiliary routines
- `fft1.f90`     → 1D FFT implementation (for spectral derivatives)


## Physical Problem

The viscous 1D Burgers equation is

$$
\frac{\partial u}{\partial t}
+ u\,\frac{\partial u}{\partial x}
= \nu\,\frac{\partial^2 u}{\partial x^2}
$$

where:

- $u(x,t)$ is the velocity field
- $\nu$ is the kinematic viscosity

This equation is a **minimal model of turbulence**:

- the nonlinear term $u\,\partial_x u$ produces **steepening** and shock formation
- the viscous term $\nu\,\partial_x^2 u$ smooths singularities

In the limit $\nu \to 0$, sharp shocks appear; for $\nu > 0$, smooth boundary
layers are obtained.

The problem is solved on a **periodic domain**:

$$
x \in [0, L]
$$

with

$$
u(x+L,t) = u(x,t).
$$

This choice is not aesthetic: it allows the use of **spectral methods (FFT)**
to compute spatial derivatives with very high accuracy.

---

## Spatial Discretization: Spectral Approach

The function $u(x,t)$ is expanded in Fourier modes:

$$
u(x,t) = \sum_k \hat{u}_k(t)\, e^{ikx}.
$$

the spatial differentiation is trivial in Fourier space

$$
\partial_x u \;\leftrightarrow\; ik\,\hat{u}_k,
$$

$$
\partial_x^2 u \;\leftrightarrow\; -k^2\,\hat{u}_k.
$$

1. FFT: $u(x) \rightarrow \hat{u}(k)$
2. Apply derivatives by multiplying by $ik$ or $-k^2$
3. IFFT: return to real space

This is what `fft1.f90` implements.

---

##  Treatment of the Nonlinear Term

The nonlinear term is evaluated in a **pseudo-spectral** way:

1. compute $u(x)$
2. compute $\partial_x u(x)$
3. multiply pointwise: $u\,\partial_x u$

Mathematically, this corresponds to a convolution in Fourier space, but
performing it in real space is far more efficient.

Note: this method can introduce **aliasing**. For moderate viscosities this is
not critical; for fully developed turbulence studies, *dealiasing* (2/3 rule)
is usually applied.

---

##  Time Discretization

The equation is reduced to a system of ODEs:

$$
\frac{d u}{dt}
= -u\,\partial_x u + \nu\,\partial_x^2 u,
$$

which is integrated in time using an explicit scheme (e.g. Euler or Runge–Kutta,
as defined in `burgers.f90`).

Conceptually:

$$
u^{n+1} = u^n + \Delta t\, \text{RHS}(u^n).
$$

Stability is controlled by:

- the CFL condition (nonlinear term)
- viscous diffusion (term $\nu k^2$)

---

## Compilation and Execution

Example using `gfortran`:

bash

- gfortran -c mod.f90 fft1.f90 burgers.f90 
- gfortran mod.f90 fft1.f90 burgers.f90

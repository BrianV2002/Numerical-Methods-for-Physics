# Spectral and Pseudo-Spectral Methods: 1D Poisson Equation

This subdirectory contains numerical implementations for solving the **one-dimensional Poisson equation with periodic boundary conditions**, comparing:

- Standard finite differences  
- A compact **Padé-type scheme** (pseudo-spectral)

The main goal is to illustrate the connection between spectral methods and compact finite-difference schemes

---

### Physical Problem

We solve the equation

\[
\frac{d^2 \phi(x)}{dx^2} = \rho(x), \qquad x \in [0,L], \qquad \phi(x+L)=\phi(x)
\]

with a known periodic source:

\[
\rho(x) = - (m k_0)^2 \sin(m k_0 x), \qquad k_0 = \frac{2\pi}{L}
\]

The exact solution is

\[
\phi_{\text{exact}}(x) = \sin(m k_0 x)
\]

which allows a direct and unambiguous measurement of numerical error.

---

### General Numerical Idea

Under periodic boundary conditions, the second-derivative operator has a **spectral interpretation**:

$$
\frac{d^2}{dx^2} \;\longrightarrow\; -k^2
$$

In a fully spectral method (FFT-based), the Poisson equation becomes algebraic in Fourier space.  
In this project, instead of using FFTs explicitly, we employ **compact Padé-type schemes**, which:

- Use a local stencil  
- Accurately approximate the spectral symbol of the derivative operator  
- Achieve near-spectral accuracy for smooth functions  

This makes them a natural example of **pseudo-spectral methods**.

---

### `rutinas.f90`

Contains linear-algebra routines for **periodic tridiagonal systems**:

- `s3p_fa`: factorizes a periodic tridiagonal matrix  
- `s3p_sl`: solves the factorized system  

These routines implement a Sherman–Morrison-type decomposition, reducing the periodic problem to a nearly tridiagonal system plus a scalar correction.

Mathematically, they solve systems of the form:

$$
A \phi = b
$$

where \(A\) represents the discrete periodic Laplacian.

---

### `poisson_diff.f90`

Solves the Poisson equation using **standard central finite differences**:

\[
\phi_{i+1} - 2\phi_i + \phi_{i-1} = dx^2\,\rho_i
\]

Characteristics:

- Second-order accuracy \(O(dx^2)\)  
- Poor representation of high-frequency Fourier modes  
- Large numerical error for oscillatory solutions  

This implementation serves as a **baseline reference**.

---

### `poisson_pad.f90`

Implements a **compact Padé scheme** for the second derivative:

\[
\frac{1}{10}\phi''_{i-1} + \phi''_i + \frac{1}{10}\phi''_{i+1}
= \frac{6}{5}\frac{\phi_{i+1}-2\phi_i+\phi_{i-1}}{dx^2}
\]

This scheme:

- Is implicit  
- Leads to a periodic tridiagonal linear system  
- Accurately reproduces the spectral operator \(-k^2\)  

For this reason, it is classified as a **pseudo-spectral method**.

---

### Compilation

bash

- gfortran rutinas.f90 poisson_diff.f90 -o poisson_diff
- gfortran rutinas.f90 poisson_pad.f90 -o poisson_pad


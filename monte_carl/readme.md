# Monte Carlo Methods

### Simple Monte Carlo Integration

Monte Carlo methods provide a numerical estimate of integrals by random sampling.
For a function \( f(x) \) defined over a domain \( V \) in \( d \) dimensions,

\[
I = \int_V f(x)\, dV,
\]

the Monte Carlo estimator using \( N \) uniformly distributed points \( x_i \) is

\[
I \approx V \langle f \rangle,
\qquad
\langle f \rangle = \frac{1}{N}\sum_{i=1}^N f(x_i).
\]

The statistical error behaves as

\[
\Delta I \sim \frac{1}{\sqrt{N}},
\]

independently of the dimension.  
This makes Monte Carlo methods inefficient for low-dimensional integrals, but
extremely powerful for high-dimensional problems.

---

### Rejection Sampling

Rejection sampling is a general method to generate random variables with a target
probability density \( p(x) \), even when direct inversion is not possible.

Given a comparison function \( f(x) \ge p(x) \), samples drawn from \( f(x) \)
are accepted with probability \( p(x)/f(x) \). Accepted values follow the
desired distribution exactly.

This method is commonly used to generate nontrivial random variates such as
Gamma or Poisson distributions.

---

### Adaptive Monte Carlo: VEGAS

The VEGAS algorithm is an adaptive Monte Carlo method based on **importance
sampling**. It improves efficiency by concentrating samples in regions where
the integrand contributes most to the variance.

VEGAS assumes that the optimal sampling distribution can be approximated by a
product of one-dimensional functions,

\[
g(x_1,\dots,x_d) = \prod_{i=1}^d g_i(x_i),
\]

which allows efficient adaptation without exponential cost.  
The sampling grid is refined iteratively, and the final result is obtained as a
weighted average over independent iterations, with consistency monitored through
a reduced \(\chi^2\).

VEGAS is particularly effective for integrands aligned with coordinate axes,
but less efficient when strong correlations are present.

---

## Monte Carlo for Statistical Physics: Ising Model

Monte Carlo techniques can also be used to simulate **statistical systems**.
A paradigmatic example is the **Ising model**, where the Metropolis algorithm
is employed to generate spin configurations according to the Boltzmann
distribution.

In this approach, trial updates are accepted or rejected based on the energy
change \(\Delta E\), ensuring convergence to thermal equilibrium. This allows
the study of phase transitions, magnetization, and critical behavior in lattice
systems.



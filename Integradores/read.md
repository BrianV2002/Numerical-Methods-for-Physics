# Second Member, Euler and RK2 Methods

This repository contains numerical solvers for ordinary differential equations
(ODEs) based on the explicit evaluation of the **second member** of the system.

The methods presented here are applied to different dynamical systems 

---

### Second Member of an ODE

A system of ODEs can be written as

$$
\frac{d\mathbf{y}}{dt} = \mathbf{F}(t,\mathbf{y})
$$

The function \(\mathbf{F}(t,\mathbf{y})\) is called the **second member**.  
It defines the dynamics of the system and contains all physical or biological
interactions (forces, damping, nonlinearities, population coupling, etc.).

Numerical integration consists of evaluating the second member in order to
advance the solution in time.

---

### Euler Method

The Euler method advances the solution using a single evaluation of the second
member at the beginning of the time step:

$$
\mathbf{y}_{n+1} = \mathbf{y}_n + \Delta t \, \mathbf{F}(t_n,\mathbf{y}_n)
$$

- First-order accurate
- Assumes the derivative remains constant over the time step
- Simple but inaccurate for nonlinear or oscillatory systems

---

### RK2 Method (Midpoint)

The second-order Runge–Kutta method improves Euler by evaluating the second
member at an intermediate state:

$$
\begin{aligned}
k_1 &= \mathbf{F}(t_n,\mathbf{y}_n) \\
k_2 &= \mathbf{F}\!\left(t_n+\tfrac{\Delta t}{2},
                         \mathbf{y}_n+\tfrac{\Delta t}{2}k_1\right) \\
\mathbf{y}_{n+1} &= \mathbf{y}_n + \Delta t\, k_2
\end{aligned}
$$

- Second-order accurate
- Uses the slope at the midpoint of the interval
- Provides better stability and phase accuracy

---

## Applications

These methods are used to integrate:
- The simple pendulum
- The forced and damped pendulum
- Predator–prey population dynamics

In all cases, the improvement from Euler to RK2 comes from a more accurate
sampling of the second member within each time step.



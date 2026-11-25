# Model conventions

## Degrees of freedom and Hamiltonian

Each site retains the six states of a dimensionless angular momentum
$J=5/2$, in ascending $J_z$ order. The four-site Hilbert-space dimension is
$6^4=1296$. Dipolar exchanges are antiferromagnetic for positive coefficients.

With each physical bond counted once, the intended lattice Hamiltonian is

$$
\begin{aligned}
H={}&J_0\sum_{\langle ij\rangle_{ab}}
  (J_i^xJ_j^x+J_i^yJ_j^y+\Delta J_i^zJ_j^z)
 +J_1\sum_{\langle ij\rangle_c}\mathbf J_i\cdot\mathbf J_j
 +J_2\sum_{\langle\langle ij\rangle\rangle_c}\mathbf J_i\cdot\mathbf J_j\\
&+\sum_{n=0}^2 K_n\sum_{\langle ij\rangle_n}
  (Q_{x,i}Q_{x,j}+Q_{y,i}Q_{y,j})\\
&+\sum_i(B_{20}O_{20,i}+B_{40}O_{40,i}+B_{44}O_{44,i})
 -g\mu_B\sum_i\mathbf B\cdot\mathbf J_i .
\end{aligned}
\tag{1}
$$

Here the quadrupolar bond classes $n=0,1,2$ are respectively nearest
in-plane, nearest-plane, and next-nearest-plane neighbors. The input uses
lowercase names such as `j0`, `k1`, and `b20`. The quadrupoles are

$$
Q_x=\frac{J_xJ_z+J_zJ_x}{2},\qquad
Q_y=\frac{J_yJ_z+J_zJ_y}{2}.
\tag{2}
$$

The Stevens-operator convention is

$$
\begin{aligned}
O_{20}&=3J_z^2-J(J+1),\\
O_{40}&=35J_z^4-[30J(J+1)-25]J_z^2
       +3J^2(J+1)^2-6J(J+1),\\
O_{44}&=\tfrac12(J_+^4+J_-^4).
\end{aligned}
\tag{3}
$$

The crystal-field coefficients therefore multiply these operators directly,
without extra Stevens factors or a change of normalization.

## Cluster and environment

The four exact plaquette bonds are `(1,2)`, `(1,3)`, `(2,4)`, and `(3,4)`.
Sites 1 and 4 belong to sublattice A; sites 2 and 3 to B.
Each site has two additional in-plane nearest neighbors represented by
self-consistent fields. The first and second neighboring planes each
contribute one neighbor above and one below on the same sublattice.

Order parameters repeat after `Nc` planes. Periodic indexing of the order
parameters preserves the physical coordination even for `Nc=1` or `Nc=2`:
coincident labels represent periodically repeated neighbors, not the removal
of physical bonds. Up to eight periods are compared, with `1 <= Nc <= 80`.

Let $\mathbf m_{a,p}$ and $\mathbf q_{a,p}$ be the dimensionless dipolar and
two-component quadrupolar expectations for sublattice $a$ in plane $p$.
For the opposite sublattice $\bar a$, the exchange fields are

$$
\begin{aligned}
\mathbf h_{a,p}
 &=2J_0(m^x_{\bar a,p},m^y_{\bar a,p},\Delta m^z_{\bar a,p})
   +J_1(\mathbf m_{a,p-1}+\mathbf m_{a,p+1})
   +J_2(\mathbf m_{a,p-2}+\mathbf m_{a,p+2}),\\
\mathbf h^{Q}_{a,p}
 &=2K_0\mathbf q_{\bar a,p}
   +K_1(\mathbf q_{a,p-1}+\mathbf q_{a,p+1})
   +K_2(\mathbf q_{a,p-2}+\mathbf q_{a,p+2}).
\end{aligned}
\tag{4}
$$

The dipolar coefficient in the cluster Hamiltonian is
$\mathbf h_{a,p}-g\mu_B\mathbf B$; the quadrupolar coefficient is
$\mathbf h^{Q}_{a,p}$. The sublattice order parameters are averages over the
two sites of each sublattice in a plaquette.

At a converged fixed point, the reported free energy per spin is

$$
f=\frac{1}{4N_c}
\left[
 \sum_{p=1}^{N_c}F^{\rm cluster}_p
 -\sum_{p=1}^{N_c}\sum_{a=A,B}
  \left(\mathbf h_{a,p}\cdot\mathbf m_{a,p}
       +\mathbf h^{Q}_{a,p}\cdot\mathbf q_{a,p}\right)
\right].
\tag{5}
$$

The subtraction removes mean-field double counting; the factor of two
sites per sublattice has already been included. At zero temperature,
$F^{\rm cluster}_p$ is the cluster ground-state energy. At positive temperature
it is obtained from the cluster partition function.

## Units and field geometry

All energy coefficients are in meV, temperature in K, and field in T.
The embedded constants are `g = 6/7`,
`mu_B = 0.05788381806 meV/T`, and `k_B = 0.08617333262 meV/K`.
The angle `theta_deg` is measured from positive $z$:

$$
\mathbf B=
\begin{cases}
B(\sin\theta,0,\cos\theta), & \texttt{inplane='x'},\\
B(\sin\theta/\sqrt2,\sin\theta/\sqrt2,\cos\theta),
 & \texttt{inplane='diag'}.
\end{cases}
\tag{6}
$$

`*_moments.csv` contains dimensionless angular-momentum expectations.
The summary moment, in Bohr magnetons per spin, is

$$
\mathbf M/\mu_B=
\frac{g}{2N_c}\sum_{p=1}^{N_c}
(\mathbf m_{A,p}+\mathbf m_{B,p}).
\tag{7}
$$

## Fourier observables

The reported dipolar wavevector is selected from the Fourier amplitude of
the in-plane staggered sublattice combination
$(\mathbf m_{A,p}-\mathbf m_{B,p})/2$. Here “in-plane staggered” refers to
the A/B alternation; all three spin components enter its amplitude.

The quadrupolar wavevector uses the sublattice-uniform combination
$(\mathbf q_{A,p}+\mathbf q_{B,p})/2$. Both searches consider discrete
stacking wavevectors $k/N_c$ from zero to one half, in reciprocal lattice
units. These observables do not survey every symmetry channel. A wavevector
attached to a numerically vanishing amplitude has no physical significance.

## Numerical scope

At zero temperature the code uses restarted Lanczos diagonalization with
full reorthogonalization. At positive temperature it uses dense Hermitian
diagonalization through LAPACK. Damped and Anderson-accelerated iterations
are available for the self-consistency loop. The program can continue
distinct branches, introduce deterministic perturbed starts, and perform
both field-sweep directions.

The residual is the maximum absolute change in the dipolar and quadrupolar
order parameters in the fixed-point update. Convergence is conditional on
the selected ansatz and numerical tolerances. It does not certify global
minimization across all periods or configurations. This package implements
the magnetic and quadrupolar solver; it does not contain an ultrasound
response calculation.

# YaHMC simulation code
**YaHMC** is a code for numerical simulations of SU(N) gauge theories with fermions in different representations. The main features of the code include:

* 4-dimensional parallelization with MPI and/or OpenMP
* Wilson-Clover fermions and Symanzik gauge action
* Minimal norm integrators (2nd and 4th order)
* Even-odd and mass preconditioned actions
* Rational approximations for single flavours
* Supported inverters
  * BiCGstab
  * CG with multishift support
  * MINRES
* Supported fermion representations
  * Fundamental
  * Two-index symmetric
  * Two-index antisymmetric
  * Adjoint
* Supported boundary conditions
  * Periodic
  * Open
  * SF
* Default observables
  * Wilson flow
  * Meson two-point functions
  * Polyakov loops

## Usage
The code has several compile flags in *programs/Makefile* that determine the simulated model. In this file it is possible to specify the number of colors, the fermion representation, the use of even/odd preconditioning, the use of clover improvement, the type of boundary conditions, the use of MPI parallelization, and the use of OpenMP. The main simulation program accepts two options:

```sh
./hmc -i config_file -o log_file
```

When running the program without any options, the configuration will be read from *hmc.cfg* and the log is written to *stdout*. The table below describes the possible settings in the configuration file, but in some cases additional knowledge of the code might be necessary to fully understand the consequences of the variables, especially in relation to the integration scheme and the action.

| Section       | Variable     | Description |
| ------------- | ------------ | ----------- |
| lattice       | dim_t        | Global dimension of lattice in the T direction. |
| lattice       | dim_x        | Global dimension of lattice in the X direction. |
| lattice       | dim_y        | Global dimension of lattice in the Y direction. |
| lattice       | dim_z        | Global dimension of lattice in the Z direction. |
| parallel      | threads      | Number of threads used for OpenMP parallelization. A value of 0 means that the default number of threads will be used. Ignored when OpenMP is disabled. |
| parallel      | np_t         | Number of MPI processes in the T direction. Ignored when MPI is disabled. |
| parallel      | np_x         | Number of MPI processes in the X direction. Ignored when MPI is disabled. |
| parallel      | np_y         | Number of MPI processes in the Y direction. Ignored when MPI is disabled. |
| parallel      | np_z         | Number of MPI processes in the Z direction. Ignored when MPI is disabled. |
| log           | level        | Specifies the amount of information written to the log file (higher number means more information). |
| rand          | seed         | Seed for the random generator. A value of 0 means that a random seed is chosen automatically. |
| cnfg          | dir          | Directory where the configurations should be saved. |
| cnfg          | prefix       | Prefix used for the configuration file names. |
| cnfg          | freq         | Specifies how often the configurations should be saved. A value of 0 means never. |
| cnfg          | start        | Starting configuration for the gauge field. Either `random` or a `unit` or the filename of an existing configuration to load. |
| cnfg          | last         | Number of configurations generated before the program exits. |
| clover        | csw          | Improvement coefficient for the clover term. Ignored when clover term is disabled. |
| clover        | cf           | Improvement coefficient used in the Dirac operator for Open and SF boundary conditions. Only used when clover improvement is enabled. |
| traj          | length       | Trajectory length for each configuration. |
| observables   | plaq_freq    | Specifies how often the Plaquette should be measured. A value of 0 means never. |
| observables   | poly_freq    | Specifies how often the Polyakov loops should be measured. A value of 0 means never. |
| observables   | mes_freq     | Specifies how often the meson two-point functions should be measured. A value of 0 means never. |
| observables   | mes_hits     | Number of sources used for calculating the meson two-point functions. |
| observables   | mes_prec     | Precision used for inverting the propagator when calculating the meson two-point functions. |
| observables   | mes_mass     | Mass used for the propagator when calculating the meson two-point functions. |
| observables   | mes_method   | Use either point sources (2) or stochastic Z2 semwall sources (1) when calculating the meson two-point functions. |


### Action
The action and the integrator structure is also defined in the configurationf file. Below is an example with a two-level integration scheme where the gauge is being integrated on the inner level and the fermions on the outer level. The code is currently written such that the gauge must be integrated on level 0 and the fermions on level 1 or higher.

```
[integrator]
level = 0
type = o2mn
steps = 1

[integrator]
level = 1
type = o2mn
steps = 10

[monomial]
level = 0
type = gauge
beta = 5.6
c0 = 1.0

[monomial]
level = 1
type = hmc
mass = -0.750
dm = 0
prec = 1e-14
mre_past = 5
```

#### Integrator
The integrator section has the following variables.

| Variable     | Description |
| ------------ | ----------- |
| level        | Integrator level with 0 being the inner most level. |
| type         | Options are `o2lf` for leapfrog and `o2mn` and `o4mn` for the 2nd and 4th order minimal norm integrators. |
| steps        | Number of integration steps. |


#### Monomials
The `gauge` monomial has the following variables.

| Variable     | Description                               |
| ------------ | ----------------------------------------- |
| type         | Monomial type `gauge`.                    |
| level        | Integrator level for this monomial.       |
| beta         | Bare coupling for the gauge field.        |
| c0           | Parameter for the improved gauge action.  |

The `hmc` monomial simulates two mass degenerate fermions.

| Variable     | Description                                                  |
| ------------ | ------------------------------------------------------------ |
| type         | Monomial type `hmc`.                                         |
| level        | Integrator level for this monomial.                          |
| mass         | Bare mass for the fermion doublet.                           |
| dm           | Mass shift when combined with the `hasenbusch` monomial.     |
| prec         | Inverter precision.                                          |
| mre_past     | Number of past solutions used in the MRE algorithm.          |

The `hasenbusch` monomial is used in combination with the `hmc` monomial as a way to precondition the action.

| Variable     | Description                               |
| ------------ | ----------------------------------------- |
| type         | Monomial type `hasenbusch`.               |
| level        | Integrator level for this monomial.       |
| mass         | Bare mass for the fermion doublet.        |
| dm           | Mass shift for the `hmc` monomial.        |
| prec         | Inverter precision.                       |


The `rhmc` monomial simulates a single fermion using rational approximations.

| Variable     | Description                               |
| ------------ | ----------------------------------------- |
| type         | Monomial type `rhmc`.                     |
| level        | Integrator level for this monomial.       |
| mass         | Bare mass for the fermion.                |
| prec         | Inverter precision.                       |
| rprec        | Precision for the rational approximation. |


## Comments
The following things might be useful to know before using the simulation code.

* Mixing OpenMP and MPI does not yield the desired performance. When running simulations on a single node, the use of OpenMP works well, but when using multiple nodes, it is recommended to use MPI only, as this gives a significantly higher performance.
* The code is not optimised when using real representations (such as the Adjoint representation). Even when the representation is real, the matrices in the Dirac operator are complex (but with the imaginary part being zero).
* It is possible to run some check routines by adding the `--run-checks` argument to the main simulation program. The program exits when the checks have finished.
* All inverter precisions are specified in terms of the squared relative error, hence the maximal precision is around 1e-30.
* When using SF boundary conditions, in the configuration file, the time direction should be T+2 to accommodate the boundary time slice.


## Disclaimer
The code is provided "as-is" without warranty of any kind. Under no circumstance does the author guarantee that the code produce useful or even correct results. However, if you discover a bug in the code, feel free to report it.

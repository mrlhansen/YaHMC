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

$ ./hmc -i config_file -o log_file

When running the program without any options, the configuration will be read from *hmc.cfg* and the log is written to *stdout*. The default simulation program will perform either a two-flavour simulation or a pure gauge simulation. The table below describes the possible settings in the configuration file, but in some cases additional knowledge of the code might be necessary to fully understand the consequences of the variables, especially in relation to the integration scheme and the action.

| Variable name | Allowed values | Description |
| ------------- | -------------- | ----------- |
| lat:dim_t     | [integer]      | Global dimension of lattice in the T direction. |
| lat:dim_x     | [integer]      | Global dimension of lattice in the X direction. |
| lat:dim_y     | [integer]      | Global dimension of lattice in the Y direction. |
| lat:dim_z     | [integer]      | Global dimension of lattice in the Z direction. |
| mp:threads    | [integer]      | Number of threads used for OpenMP parallelization. A value of 0 means that the default number of threads will be used. Ignored when OpenMP is disabled. |
| mp:np_t       | [integer]      | Number of MPI processes in the T direction. Ignored when MPI is disabled. |
| mp:np_x       | [integer]      | Number of MPI processes in the X direction. Ignored when MPI is disabled. |
| mp:np_y       | [integer]      | Number of MPI processes in the Y direction. Ignored when MPI is disabled. |
| mp:np_z       | [integer]      | Number of MPI processes in the Z direction. Ignored when MPI is disabled. |
| log:level     | 0-40           | Specifies the amount of information written to the log file (higher number means more information). |
| rand:seed     | [integer]      | Seed for the random generator. A value of 0 means that a random seed is chosen automatically. |
| run:start     | random,unit,[string] | Starting configuration for the gauge field. It is possible to start from either a random or a unit gauge configuration, or to load an existing configuration by specifying the file name of the configuration here. |
| run:beta      | [float]        | Bare coupling for the gauge field. |
| run:c0        | [float]        | Parameter for the gauge action. |
| run:mass      | [float]        | Bare mass for the fermion doublet. Ignored for pure gauge simulations. |
| run:csw       | [float]        | Improvement coefficient for the clover term. Ignored when clover term is disabled and for pure gauge simulations. |
| run:last      | [integer]      | Number of configurations generated before the program exits. |
| act:puregauge | 0,1            | A non-zero value indicates that only the gauge field should be simulated. |
| act:dm        | [float]        | Mass shift used in the  mass preconditioned action. When the value is zero, mass preconditioning is disabled. |
| save:dir      | [string]       | Directory where the configurations should be saved. |
| save:name     | [string]       | Prefix used for the configuration file names. |
| save:freq     | [integer]      | Specifies how often the configurations should be saved. A value of 0 means never. |
| plaq:freq     | [integer]      | Specifies how often the Plaquette should be measured. A value of 0 means never. |
| poly:freq     | [integer]      | Specifies how often the Polyakov loops should be measured. A value of 0 means never. |
| mes:freq      | [integer]      | Specifies how often the meson two-point functions should be measured. A value of 0 means never. |
| mes:hits      | [integer]      | Number of sources used for calculating the meson two-point functions. |
| mes:prec      | [float]        | Precision used for inverting the propagator when calculating the meson two-point functions. |
| mes:method    | 1,2            | Use either point sources (2) or stochastic Z2 semwall sources (1) when calculating the meson two-point functions. |
| int:nsteps    | [integer]      | Number of integration steps used on the coarse fermion level. Ignored for pure gauge simulations. |
| int:hsteps    | [integer]      | Number of integration steps used on the fine fermion level. Only used for the mass preconditioned action. |
| int:gsteps    | [integer]      | Number of integration steps used for the gauge field. |
| int:length    | [float]        | Integration length for each configuration. |
| inv:prec      | [float]        | Precision used for inverting the Dirac operator during the trajectory. |
| mre:past      | [integer]      | Number of past solutions used in the chronological solver. |
| bc:cf         | [float]        | Improvement coefficient used in the Dirac operator for Open and SF boundary conditions. Only used when clover improvement is enabled. |

## Comments
The following things might be useful to know before using the simulation code.

* Mixing OpenMP and MPI does not yield the desired performance. When running simulations on a single node, the use of OpenMP works well, but when using multiple nodes, it is recommended to use MPI only, as this gives a significantly higher performance.
* The default action and integrator structure is programmed in the *library/update.cpp* file. Have a look here for details about the action.
* The code is not optimised when using real representations (such as the Adjoint representation). Even when the representation is real, the matrices in the Dirac operator are complex (but with the imaginary part being zero).
* It is possible to run some check routines by adding the "--run-checks" argument to the main simulation program. The program exits when the checks have finished.
* All inverter precisions are specified in terms of the squared relative error, hence the maximal precision is around 1e-30.
* When using SF boundary conditions, in the configuration file, the time direction should be T+2 to accommodate the boundary time slice.
* At the moment the code is written such that the gauge must be integrated on level 0 and the fermions must be integrated on level 1 or higher.

## Disclaimer
The code is provided "as-is" without warranty of any kind. Under no circumstance does the author guarantee that the code produce useful or even correct results. However, if you discover a bug in the code, feel free to report it.

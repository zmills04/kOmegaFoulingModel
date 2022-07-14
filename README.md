# kOmegaFoulingModel
This code simulates the development of a fouling layer along a single
channel of an exhaust gas recirculation heat exchanger. It utilizes the
Lattice Boltzmann Method (LBM) for modeling fluid flow, a finite differences
for modeling the heat transport, and Brownian Dynamics (BD) for modeling the
motion of sub-micron sized soot particles. A sticking probably based boundary
condition allows the BD model to incorporate the effects of shear along the wall
with particles being more likely to deposit in locations with low shear, and less
likely for higher shear. The deposited particles in turn cause the fouling layer to
grow, which is modeled by shifting wall surfaces into the fluid. A detailed description
of the implementation can be found in
[my thesis](https://smartech.gatech.edu/bitstream/handle/1853/58218/MILLS-DISSERTATION-2017.pdf).

The majority of the methods in this code has been implemented in OpenCL to allow for GPU
acceleration, which the Lattice Boltzmann Method is well suited for. Furthermore, a sorting
algorithm and multiple linear solvers obtained from the AMD Bolt and clSparse libraries,
respectively, have been directly implemented in the code as well. This was done as a
learning exercise and to avoid unnecessary CPU-GPU data transfers since
all operations were performed on the GPU and no transfers were necessary.

The original code used for my thesis work was only capable of simulating laminar flow.
Several attempts were made to implement the capability to model turbulent flows. Initially
several Entropic Lattice Boltzmann Methods were implemented in the code. The first ELBM method
implemented used a single entropic relaxation parameter, alpha, to enforce the second law of
thermodynamics. With this implemented, the model was capable of simulating flows with Re > 2000,
but the model did not prouduce a turbulent velocity profile in a straight channel in a 2D simulation
(at the time I did not have access to HPC resources, so I was limited to a single 8GB GPU, not nearly
enough for a large 3D simulation). Assuming that it might be a limitation of the single relaxation parameter
model, I implemented the multi-relaxation time ELBM developed by Karlin et al. This method again was stable,
but unable to capture the turbulent velocity profile in a straight channel. After discussing the issues with
several LBM experts, I learned that the ELBM cannot produce a turbulent velocity profile in 2 dimensions.
After abandoning the ELBM method, I decided to implement a k-omega turbulence model in the code.
I was successful in getting it to work in a straight channel and a wavy channel. Unfortunately,
the model was not stable when the effects of the growing fouling layer were included in the simulations.

At this point I have decided that the lattice Boltzmann method is not the right tool for this system and have not
continued beyone the attempt to implement a k-omega turbulence model. I have uploaded the code so that anyone
interested in implementing the lattice Boltzmann method in OpenCL can see an example. All openCL kernels can be
found in source/Kernels. There is also an implementation of the ELBM and an implementation of LBM which uses a few tricks
to avoid loss of accuracy when using floats instead of doubles to track the distributions. These can be found in 
source/Kernels/oldKernels.

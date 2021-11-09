# Lattice Boltzmann Fluid Simulation

The Lattice Boltzmann Method (LBM) is a numerical method for fluid simulation where, instead of solving the Navier-Stokes equations directly, the motion of the fluid is recovered from the velocity distribution of the constituting molecules and approximations of molecular collisions.
Simulations using LBM are implemented on a discretized domain accompanied by a discretized velocity distribution, with a finite number of possible velocity magnitudes and directions. The updating of the velocity distributions can be separated in two distinct stages:

* Collision step. The velocity distribution at each individual cell is adjusted towards the equi- librium distribution, simulating collisions between molecules.

* Streaming step. Following the collision step, the velocity distribution at neighboring cells are updated based on the velocity distribution at the current cell, simulating movement of molecules in space.

Since the dynamical rules in the method is local (each cell is only updated in relation to neighboring cells) it can be made highly efficient on parallel computers.

This repository contains a parallel implementation of the multiphase multiphase Shan-Chen Lattice-Boltzmann method suitable for massively parallel computing, using the MPI standard.

![animated example](https://github.com/Haspaker/parallel-lattice-boltzmann-simulation/blob/main/example.gif?raw=true)

The above example shows multiphase vapor-liquid coexistence, with droplets coalescing from small initial fluctuations around a mean value.
# AFAC_solver
Asynchronous Fast adaptive composite Grid solver for partial differential equations allowing variational formulation. Writen in Fortran using Coarrays for distributed memory parallelism. 

## Implementation

This is a proof of concept code it shows the working of the AFAC algorithm on a nested grid with an arbitrary number of nested levels. It calculates the gravitational potential of a given density distribution and compares them to the analytical results. 
The results are outputted as .dat files. 

## Results

At the moment the following problems have working analytical solutions: 
Sphere
Condensed Sphere
Binary Spheres 

The exact output can be examined in the src/output.f90 file. 


## Contributing

Contributions are welcome! Please open an issue or submit a pull request.

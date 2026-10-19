# AFAC_solver
Asynchronous Fast adaptive composite Grid solver for partial differential equations allowing variational formulation. Writen in Fortran using Coarrays for distributed memory parallelism. 

## Implementation

This is a proof of concept code it shows the working of the AFAC algorithm on a nested grid with an arbitrary number of nested levels. It calculates the gravitational potential of a given density distribution and compares them to the analytical results. 
The results are outputted as .dat files. The amount of nested levels can be changed by changing the config.caf file. 

## Results

At the moment the following problems have working analytical solutions: 
Sphere
Condensed Sphere
Binary Spheres 

The exact output can be examined in the src/output.f90 file. 

Plots of the results will be following soon. 

The Benchmark results of a single image run can be found in benchmark_results.csv

I will update the Benchmarks soon for multi-image runs. 

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.

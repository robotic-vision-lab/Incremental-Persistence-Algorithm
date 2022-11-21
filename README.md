## Incremental Persistence Algorithm

This software implements the algorithms for computing topological persistence 
as described in [1]. Given a 3D point cloud, the software computes the zeroth 
(Betti 0) and first (Betti 1) homology groups which correspond to the number of 
connected components and holes in the dataset, respectively. This is done by first 
approximating the input space using a Vietoris-Rips complex. Next, the simplicial 
persistence pairs for a given resolution are computed incrementally. Finally, the 
birth and death times of the 0- and 1-simplices are recorded to an output file. 

### Running MATLAB

To start a MATLAB session, invoke the run\_matlab.sh script which adds the 
necessary paths and environment variables to the current session.  You may have 
to customize this script for your environment.

### How to Compute Topological Persistence 

The software takes as input a PCD file along with a maximum scale parameter.  To 
compute the Betti numbers for the dataset run the following command:

    >> ipa('in_file.pcd', 'out_file.txt', 0.003)

The birth and death times of the simplices are written to 'out\_file.txt.'

### References

[1] H. Edelsbrunner, D. Letscher, and A. Zomorodian, "Topological persistence 
    and simplification," Discrete and Computational Geometry, pp. 511-533, 2002.  
    
### License

[![License](https://img.shields.io/badge/License-BSD_2--Clause-orange.svg)](https://github.com/robotic-vision-lab/Incremental-Persistence-Algorithms/blob/master/LICENSE)

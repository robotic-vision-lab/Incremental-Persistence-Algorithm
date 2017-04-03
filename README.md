#### Incremental persistence algorithm 

This software implements the algorithms for computing topological persistence 
as described in [1].  Given a 3D point cloud as input, the software computes the
zeroth (Betti 0) and first (Betti 1) homology groups which correspond to the 
number of connected components and holes in the dataset, respectively.  This is
done by first approximating the input space using a Vietoris-Rips complex.  
Then, the simplicial persistence pairs for a given resolution are computed 
incrementally.  Finally, the birth and death times of the 0- and 1-simplices are 
recorded to an output file. 

#### Running MATLAB

To start a MATLAB session, invoke the run\_matlab.sh script which adds the 
necessary paths and environment variables to the current session.  You may have 
to customize this script for your environment.

#### How to compute topological persistence 

The software takes as input a PCD file along with a maximum scale parameter.  To 
compute the Betti numbers for the dataset run the following command:

    >> ipa('in_file.pcd', 'out_file.txt', 0.003)

The birth and death times of the simplices are written to 'out\_file.txt'.

#### Notes

In matpcl, there's a bug in the splitting of the RGB values.  The shift needs to
be to the right:

    $ diff loadpcd.m.orig loadpcd.m
    231,232c231,232
    <                 R = double(bitand(255, bitshift(rgb, 16))) /255;
    <                 G = double(bitand(255, bitshift(rgb, 8))) /255;
    ---
    >                 R = double(bitand(255, bitshift(rgb, -16))) /255;
    >                 G = double(bitand(255, bitshift(rgb, -8))) /255;
    237,239c237,239
    <                 R = double(bitand(255, bitshift(rgb, 24))) /255;
    <                 G = double(bitand(255, bitshift(rgb, 16))) /255;
    <                 B = double(bitand(255, bitshift(rgb, 8))) /255;
    ---
    >                 R = double(bitand(255, bitshift(rgb, -24))) /255;
    >                 G = double(bitand(255, bitshift(rgb, -16))) /255;
    >                 B = double(bitand(255, bitshift(rgb, -8))) /255;

#### References

[1] H. Edelsbrunner, D. Letscher, and A. Zomorodian, "Topological persistence 
    and simplification", Discrete and Computational Geometry, pp. 511-533, 2002.    

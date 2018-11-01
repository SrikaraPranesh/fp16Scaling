# Scaling into fp16
MATLAB codes for scaling a matrix given in fp64 or fp32
into fp16, and using the scaled matrix to solve system
of linear equations using GMRES based iterative refinement.

## Related publications
* N. J. Higham, S. Pranesh, M. Zounon. [Squeezing a Matrix Into Half Precision, with 
an Application to Solving Linear Systems]

## Main Execution file
* **_MAIN.m_** Generates all the tables in the manuscript into the folder results/.


## Requirements
* The codes have been developed and tested with MATLAB 2014b.
* This code requires Cleve Laboratory to perform half precision computations and 
Advanpix Multiprecision Computing Toolbox for extended precision computations. 
A free trial of Advanpix is available for download from https://www.advanpix.com/.
  

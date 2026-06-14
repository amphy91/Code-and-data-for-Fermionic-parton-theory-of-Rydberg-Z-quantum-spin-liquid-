# Code-and-data-for-Fermionic-parton-theory-of-Rydberg-Z-quantum-spin-liquid-
This repository contains the codes and data files for the paper "Fermionic parton theory of Rydberg Z₂ quantum spin liquid"

Contents:
1. The file "Ansatze_Table-I.nb" contains the mathematica notebook for calculating the symmetry allowed mean-field parameters given in Table I using the PSG solutions given in the supplimental information.
2. The file "dsf.f90" contains the code for calculating the data for Fig.3. The corresponding data file is given by "dsf.dat".
3. The file "ssf.f90" contains the code for calculating the data for Fig.2(b). The corresponding data file is given by "ssf.dat".
4.There are two mathematica notebooks "DSF.nb" and "SSF.nb" which produce the corresponding plot with appropriate colorscheme.

Requirements:
The codes were developed using Gfortran and mathematica
Required packages:
blas and lapack

Reproducing the results:
1. Use "Ansatze_Table-I.nb" for reproducing the symmetry allowed mean-field parameters given in Table I.
2. Run "dsf.f90" and "ssf.f90" to generate the data file using the library blas and lapack.
3. Use the mathematica notebooks for plotting.

Contact:
For question regarding the code, please contact atanu.maity@uni-wuerzburg.de/amphy91@gmail.com


# PWE-JDQR
GNU-Octave / MATLAB implementation of the application of QR-type Jacobi-Davidson algorithms for interior eigenvalue calculations in plane wave expansion in photonics.

Coded by Thomas Kamalakis, Department of Informatics and Telematics, Harokopio University Greece

Tested in GNU/Octave version 4.0.0 and MATLAB 8.5.0.197613 (R2015a)

Please refer to LICENCE.txt for details on licencing.

m-file description:
  
  -PWE_JDQR.m is an example m-file that illustrates the use of the JDQR algorithm.
  -jdqr.m is my implementation of the JDQR algorithm for HERMITIAN matrices. Do not use in non-Hermitian eigen-problems.
  -MV.m describes the action of the eigen-matrix on a vector
  -dielectric_tensors.m calculates the dielectric tensors required for averaging
  -lattice_vectors_rect.m calculates the reciprocal lattice vectors.
  -dielectric_averaging2.m performs dielectric averaging.
  -mymgs.m performs modified Gram-Schmidt orthogonalization.
  -orth_vecs.m calculates the orthogonal eigenvectors used in the magnetic field expansion.
  

# project-code-git
Two fortran codes to generate a hermitian matrix, diagonlize the same and use its eigenvectors to change the basis of another non-hermitian matrix

The first code uses LAPACK library for diagonalization

The second code is a parallized code that uses MPI with block size of 32x32 and grid size of 2x2 & uses SCALAPACK for diagonalization. Note that the data files are not uploaded due to large size

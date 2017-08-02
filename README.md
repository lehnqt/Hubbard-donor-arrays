Matlab code for exact diagonalisation of an extended Hubbard model and computing the conductance of a donor array in silicon. This code was used in arXiv:1707.01876 (An extended Hubbard model for mesoscopic transport in donor arrays in silicon).

Details of the subfolders:

1. EigsMatlab: use the Matlab built-in eigs function for diagonalisation (suitable for arrays with less than 12 sites).

2. LanczosGPU: use Lanczos algorithm run on a GPU for diagonalisation (for arrays with up to 16 sites).

3. HartreeFock: The ground state at half-filling is obtained by the Hartree-Fock mean field calculation (for large two dimensional arrays).

4. HoleHubbard: The ground state near half-filling is obtained by the hole-hubbard approximation described in the Appendix of the paper (for large two dimensional arrays).

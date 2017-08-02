Matlab code for exact diagonalisation of an extended Hubbard model and computing the conductance of a donor array in silicon. This code was used in arXiv:1707.01876 (An extended Hubbard model for mesoscopic transport in donor arrays in silicon).
Details of the subfolders:
EigsMatlab: use the Matlab built-in eigs function for diagonalisation (suitable for arrays with less than 12 sites).
LanczosGPU: use Lanczos algorithm run on a GPU for diagonalisation (for arrays with up to 16 sites).

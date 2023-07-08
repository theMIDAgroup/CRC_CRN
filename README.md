# CRC_CRN

Matlab code for the in-silico simulation of the effects of mutations and targeted therapies in a chemical reaction network modeling the G1-S transition phase of a colorectal cancer cell.

We also provide the implementation of the NLPC (non-linearly projected combined) method, a fast and globally convergent combined Newton and gradient descent method for computing equilibrium points of chemical reaction networks.

The code allows to reproduce the results shown in 

* Sommariva, S., Caviglia, G., Ravera, S., Frassoni, F., Benvenuto, F., Tortolina, L., Castagnino, N., Parodi, S., Piana, M. (2021) Computational quantification of global effects induced by mutations and drugs in signaling networks of colorectal cancer cells. Scientific reports 11(1), 1-13

* Berra, S., La Torraca, A., Benvenuto, F., Sommariva, S., A fast and convergent combined Newton and gradient descent method for computing steady states of chemical reaction networks, submitted.

* Sommariva, S., Berra, S., Biddau, G., Caviglia, G., Benvenuto, F., Piana, M., In-silico modelling of the mitogen-activated protein kinase (MAPK) pathway in colorectal cancer: mutations and targeted therapy, submitted.

# For reproducing results from Sommariva et al. 2021:

Code is written in Matlab R2015b.
Results of the simulations are provided in the folder './results'.

Run script_fig_\*.m to reproduce the figures of the paper. 
Run main_script_\*.m to reproduce the results of the paper.

# For reproducing results from Berra et al. submitted:

Code is written in Matlab R2017b.
Results of the simulations are provided in the folder './test_NLPC/results_paper'.
The implementation of the NLPC methods can be found in './funcs/f_NLPC_restart.m'

To reproduce the figures of the paper, run ./test_NLPC/script_figure_\*.m .
To reproduce the analysis of the paper, run:
* ./test_NLPC/main_extract_x0.m to sample the initial points for NLPC.
* ./test_NLPC/main_NLPC.m and ./test_NLPC/main_NLPC_orthogonal.m to run NLPC with the proposed non linear projector and the classical orthogonal projector, respectively
* ./test_NLPC/main_dyn.m to compute the steady states with a classical ODEs solver.

# For reproducing results from Sommariva et al. submitted:

Code is written in Matlab R2017b.
Results of the simulations are provided in the folder './test_MAPK/results_paper'.

To reproduce the figures of the paper, run ./test_MAPK/script_figure_\*.m .
To reproduce the analysis of the paper, run the m-files ./test_MAPK/main_\*.m .The scope of each script is extensively described therein.

# External packages:

Our code makes use of the following external graphical tools:
* Robert (2021). symlog (https://github.com/raaperrotta/symlog), GitHub. Retrieved January 15, 2021.
* Povilas Karvelis (2022). daboxplot (https://github.com/frank-pk/DataViz), GitHub. Retrieved December 28, 2022.

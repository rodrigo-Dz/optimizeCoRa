using DifferentialEquations
using ForwardDiff
using LinearAlgebra
using Plots
using NLsolve
using Sobol
using DelimitedFiles;
using Statistics

# Main arguments
iARG = (mm = "ATFv1_p02",		# Label for motif file
        ex = "Fig1",		# Label for parameters file
        pp = :mY,			# Label for perturbation type
        ax = :mY,               # Label for condition/environment
        an = "Optimize"		# Label for analysis type (Explore, Optimize)
)

# Main functions
fn = include(string("Library/FN_CoRa.jl"));

# Perturbation details
pert = (p   = iARG.pp,	# Parameter to be perturbed
        d   = 1.05,		# Perturbation size (Delta rho)
        c   = iARG.ax,	# Condition parameter
        r   = [-3.0,3.0],	# Range of conditions
        eps = 0.1,
        coras = 50,             # Number of CoRa simulations
        solver = "fast",        # Solover to use ("fast", or "slow" (more precise))
        tspan = (0.0, 1e8));        

# Edit the analysis you want:

# exploration details
expl  = (n_params = 3,
        pOp  = [:mU,:mW,:eP],	# Parameters to optimize
        pMin = [-3.0, -2.0, -1.0],           # Mínimo log10 para cada uno
        pMax = [3.0,  2.0,  1.0],            # Máximo log10 para cada uno      
        n_points = 2048,           # number of points to evaluate
        prtD =1);		# flag for printing full DY curve 

# Optimization details
opt  = (n_params = 3,
        pOp  = [:mU,:mW,:eP],	# Parameters to optimize
        pMin = [-3.0, -3.0, -3.0],           # Mínimo log10 para cada uno
        pMax = [3.0,  3.0,  3.0],            # Máximo log10 para cada uno      
        iter = 10000,                            # number of points to evaluate
        cov = [0.1, 0.1, 0.1],                  # covariance for each parameter
        M = 10,                            # mutation step size
        prtD =1,
        rand = 1,
        );	


# Model
mm = include(string("Library/Md_",iARG.mm,".jl"));

# Initial Conditions
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] 

# Core parameters
include(string("InputFiles/ARGS_",iARG.mm,"_",iARG.ex,"_Par.jl"))


if (iARG.an == "Explore")
        fn.explore_all(
                iARG,       # motif     
                mm,         # model
                p,          # core parameters
                pert,       # perturbation details
                expl,        # optimization details
                u0)         # initial conditions
elseif (iARG.an == "Optimize")
        fn.optimize(
                iARG,       # motif     
                mm,         # model
                p,          # core parameters
                pert,       # perturbation details
                opt,        # optimization details
                u0)         # initial conditions
else
        println("Error: Analysis type not recognized")
end





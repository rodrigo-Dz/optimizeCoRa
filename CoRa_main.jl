using DifferentialEquations
using ForwardDiff
using LinearAlgebra
using Plots
using NLsolve
using Sobol
using DelimitedFiles;
using Statistics

# Main arguments
iARG = (mm = "ATFv1",		# Label for motif file
        ex = "Fig1",		# Label for parameters file
        pp = :mY,			# Label for perturbation type
        ax = :mY,			# Label for condition/environment
)

# Main functions
fn = include(string("Library/FN_CoRa.jl"));

# Perturbation details
pert = (p   = iARG.pp,	# Parameter to be perturbed
        d   = 1.05,		# Perturbation size (Delta rho)
        c   = iARG.ax,	# Condition parameter
        r   = [-3.0,3.0],	# Range of conditions
        eps = 0.1,
        coras = 50);        # number of CoRas for each coRa curve

# oprimization details
opt  = (n_params = 3,
        pOp  = [:mU,:mW,:eP],	# Parameters to optimize
        pMin = [-3.0, -2.0, -1.0],           # Mínimo log10 para cada uno
        pMax = [3.0,  2.0,  1.0],            # Máximo log10 para cada uno
        solver = "fast",        # Solover to use ("fast", or "slow" (more precise))
        tspan = (0.0, 1e8),      # tspan for "slow" solver
        n_points = 2048,           # number of points to evaluate
        prtD =1);		# flag for printing full DY curve  

# Model
mm = include(string("Library/Md_",iARG.mm,".jl"));

# Initial Conditions
u0 = [0.0, 0.0, 0.0, 0.0]  # Y, U, W, C


# Core parameters
include(string("InputFiles/ARGS_",iARG.mm,"_",iARG.ex,"_Par.jl"))

fn.explore_all(
        iARG,       # motif     
        mm,         # model
        p,          # core parameters
        pert,       # perturbation details
        opt,        # optimization details
        u0)         # initial conditions




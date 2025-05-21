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
        # ex = "Fig3",		# Label for parameters file
        pp = :mY,			# Label for perturbation type
        ax = :mY,			# Label for condition/environment
)
# Perturbation details
pert = (p   = iARG.pp,	# Parameter to be perturbed
        d   = 1.05,		# Perturbation size (Delta rho)
        c   = iARG.ax,	# Condition parameter
        r   = [-3.0,3.0],	# Range of conditions
        eps = 0.1,
        coras = 50);        # number of CoRas for each coRa curve

# Initial Conditions
u0 = [0.0, 0.0, 0.0, 0.0]  # Y, U, W, C

# oprimization details
opt  = (n_params = 3,
        pOp  = [:mU,:mW,:eP],	# Parameters to optimize
		pMin = -3.0,		        # Minimum parameter value to explore (log10)
		pMax = 3.0,			    # Maximum parameter value to explore (log10)       
        solver = "fast",        # Solover to use ("fast", or "slow" (more precise))
        tspan = (0.0, 1e8),      # tspan for "slow" solver
        n_points = 8);		    # number of points to evaluate

println(opt.pOp)
# Model
mm = include(string("Library/Md_",iARG.mm,".jl"));
# Main functions
fn = include(string("Library/FN_CoRa.jl"));

# Core parameters
include(string("InputFiles/ARGS_",iARG.mm,"_Par.jl"))


# range to explore
ran = 10 .^ range(-3, 3, length=pert.coras) 

holi

fn.explore_all(
    mm.ode_system!, # system with feedback
    mm.ode_systemNF!, # system without feedback
    p,      # core parameters
    pert,
    opt,
    ran,
    u0, 
    mm)
    
    



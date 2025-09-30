
# --- Main script for CoRa analysis ---
# --- Rodrigo Aguilar
# --- September 2025


# --- Required libraries ---
using DifferentialEquations
using ForwardDiff
using LinearAlgebra
using Plots
using NLsolve
using Sobol
using DelimitedFiles
using Statistics

# --- Main functions ---
fn = include(string("Library/FN_CoRa.jl"))


# --- Main arguments ---
iARG = (mm = "ATFv1_p02",       # Label for motif file
        ex = "Fig1",		# Label for parameters file
        pp = :mY,		# Label for perturbation type
        ax = :mY,               # Label for condition/environment
        an = "Explore"		# Label for analysis type (Explore, Optimize, Dynamics, Curve)
)

# --- Perturbation details ---
pert = (p   = iARG.pp,	        # Parameter to be perturbed
        d   = 1.05,		# Perturbation size (Delta rho)
        c   = iARG.ax,	        # Condition parameter
        r   = [-3.0,3.0],	# Range of conditions
        eps = 0.1,              # eps 
        coras = 30,             # Number of conditions for each curve
        solver = "fast",        # Solover to use ("fast", or "slow" (more precise))
        tspan = (0.0, 1e8)      # tspan for slow solver   
)     

# --- Model ---
mm = include(string("Library/Md_",iARG.mm,".jl"))

# --- Core parameters ---
include(string("InputFiles/ARGS_",iARG.mm,"_",iARG.ex,"_Par.jl"))




# --- Edit the analysis you want: ---
     
# Exploration details
expl  = (pOp  = [:mU,:mW,:eP],	        # Parameters to optimize
        pMin = [-3.0, -3.0, -3.0],      # Min (log scale)
        pMax = [3.0,  3.0,  3.0],       # Max (log scale) 
        n_params = length(pOp),         # number of parameters to explore    
        n_points = 4096,                # number of points to evaluate (power of 2 recommended)
        prtD =1)		        # flag for printing full curve 

# Optimization details
opt  = (n_params = 7,
        pOp  = [:mU, :mW, :eP, :gU, :gW, :e0, :eM],	        # Parameters to optimize
        pMin = [-3.0, -3.0, -3.0, -3.0, -3.0, -3.0, -3.0],      # Mínimo log10 para cada uno
        pMax = [3.0, 3.0, 3.0, 1.0, 1.0, 1.0, 1.0],             # Máximo log10 para cada uno   
        n_params = length(pOp),                                 # number of parameters to explore       
        iter = 50000,                                           # number of points to evaluate
        cov = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],              # covariance for each parameter
        M = 10,                                                 # mutation step size
        prtD =1,                                                # flag for printing full curve
        rand = 1,                                               # flag for random initial parameters
        )	


        
# --- Run ---

if (iARG.an == "Explore")
        fn.explore(
                iARG,       # Main arguments     
                mm,         # Model
                p,          # Core parameters
                pert,       # Perturbation details
                expl,       # Exploration details
                u0)         # Initial conditions
elseif (iARG.an == "Optimize")
        fn.optimize(
                iARG,       # Main arguments     
                mm,         # Model
                p,          # Core parameters
                pert,       # Perturbation details
                opt,        # Optimization details
                u0)         # Initial conditions
elseif (iARG.an == "Dynamics")
        fn.dynamics(
                mm,         # Model
                u0,         # Initial conditions
                p,          # Core parameters
                tspan,      # tspan for dynamics
                pert,       # Perturbation details
                iARG)       # initial conditions
elseif (iARG.an == "Curve")
        fn.curve(
                p,          # Core parameters
                u0,         # Initial conditions
                mm,         # Model
                pert)       # initial conditions
else
        println("Error: Analysis type not recognized")
end





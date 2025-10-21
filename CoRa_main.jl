
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

# --- Settings ---
Settings_file = "ARGS_FADv1_p01_Fig2_Settings.jl"
include(string("./InputFiles/", Settings_file))

# --- Main functions ---
fn = include(string("Library/FN_CoRa.jl"))

# --- Model ---
mm = include(string("Library/Md_",iARG.mm,".jl"))

# --- Core parameters ---
include(string("InputFiles/ARGS_",iARG.mm,"_",iARG.ex,"_Par.jl"))

# ------

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
                pert,       # Perturbation details
                dyn,        # Dynamics details
                iARG)       # Main arguments
elseif (iARG.an == "Curve")
        fn.curve(
                p,          # Core parameters
                u0,         # Initial conditions
                mm,         # Model
                pert,       # Perturbation details
                iARG)       # Main arguments
else
        println("Error: Analysis type not recognized")
end





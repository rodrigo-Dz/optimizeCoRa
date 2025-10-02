
# --- Main functions for CoRa analysis
# --- Rodrigo Aguilar
# --- September 2025


module fn
# --- Required libraries ---
	using DelimitedFiles
	using Distributions
	using DifferentialEquations
	using ModelingToolkit
	using Sobol
	using NLsolve
	using ForwardDiff
	using LinearAlgebra
    using Statistics
    using ParameterizedFunctions
    using Plots


# --- Solve for equilibrium (fast) ---
#     p ------ Set of parameters
#     u0 ----- Initial conditions
#     system - set of ODEs
#     Returns: 
#        [0] vector with equilibrium points
#        [1] error code (0= failed to converge, 1= converged)

    function find_equilibrium(p, u0, system; max_iters=1000, ftol=1e-12, xtol=1e-12)
        par = collect(values(p))
        
    # Define the equilibrium condition function   
        function equilibrium_condition(F, u)
            du = similar(u)
            system(du, u, par, 0.0)
            F .= du
            return F
        end
        
    # Primary method: Newton's method
        result = nlsolve(equilibrium_condition, u0, 
                        method=:newton,
                        autodiff=:forward,
                        ftol=ftol,
                        xtol=xtol,
                        iterations=max_iters,
                        show_trace=false)

    # If primary method fails, try trust region method with relaxed tolerances
        if !result.f_converged && !result.x_converged
            @warn "Primary method failed, trying trust region method"
            result = nlsolve(equilibrium_condition, u0,
                            method=:trust_region,
                            autodiff=:forward,
                            ftol=ftol*10,     # relaxed tolerances
                            xtol=xtol*10,
                            iterations=max_iters)
        end
        
    # Final check
        if !result.f_converged && !result.x_converged
            @warn "Equilibrium search failed to converge, returning NaN"
            return result.zero, 0     
        else
            return result.zero, 1
        end
    end



# --- Calculate Jacobian matrix at equilibrium ---
#     u_eq --- Equilibrium point
#     p ------ Set of parameters
#     system - set of ODEs
#     Returns: 
#        Jacobian matrix

    function compute_jacobian(u_eq, p, system)
        p_values = collect(values(p))
        function equilibrium_condition(u)
            du = similar(u)
            system(du, u, p_values, 0.0)
            return du
        end
        ForwardDiff.jacobian(equilibrium_condition, u_eq)
    end


# --- Evaluate CoRa ---
#     p ----- System parameters
#     u0 ---- Initial conditions
#     mm ---- Models 
#     pert -- Perturbation details
#     Returns: 
#        CoRa value, 
#        steady state of the system, 
#        steady state of the controlled variable,
#        type of error (0=no error, 1=could not find positive SS, 2=oscillations, 3=other errors)

    function evalCoRa(p, u0, mm, pert)
        copy_p = copy(p)
        error = 0

    # If initial conditions are NaN, set to zero
        if isnan(u0[1])
            u0 = fill(0.0, length(u0))
        end 

    # Solve Feedback system
        println("Solving F...")
        SS, eq = fn.find_equilibrium(p, u0, mm.FB)
        FB = mm.out_FB(SS)  

    # If negative values, try new initial conditions
        i = 1
        while any(x -> x < 0, SS)
            println("negative values in SS, trying new initial conditions")
            u_ = length(u0)
            u_ = fill(10.0^i, u_)

            SS, eq = fn.find_equilibrium(p, u_, mm.FB)
            FB = mm.out_FB(SS)  

            i += 1

            if i > 5
                println("could not find positive SS")
                error = 1
                CoRa  = NaN
                FB   = NaN
                break
            end
        end

    # If the solver failed, check for oscillations
        if eq == 0
            J = fn.compute_jacobian(SS, p, mm.FB)
            eigenvalues = eigvals(J)  
        # Complex eigenvalues with positive real part indicate oscillations
            has_oscilations = any(imag(l) != 0 && real(l) >= 0 for l in eigenvalues) 

            if has_oscilations == true
                CoRa = NaN  
                FB = NaN
                error = 2
                println("The system has oscillations")
            else
                error = 3
                CoRa  = NaN
                FB   = NaN
                println("Unknown error")
            end
        end

    # Calculate CoRa only if no errors
        if error == 0   
        # Solve NF system    
             println("Solving NF...")

            mm.localNF(p,SS)      # Adjust parameters for no feedback system
            SS_nFB, _ = fn.find_equilibrium(p, SS, mm.nFB)
            nFB = mm.out_FB(SS_nFB)
        # If relative difference between FB and nFB is less than 0.0001, proceed to perturbation   
            if (abs(FB - nFB))/ FB < 0.01  
            # perturbation to FB
                p = copy(copy_p)
                p[pert.p] = p[pert.p]*pert.d
                SS_FBp, _ = fn.find_equilibrium(p, SS, mm.FB)
                FB_p = mm.out_FB(SS_FBp)

            # perturbation to nFB    
                p = copy(copy_p)
                p[pert.p] = p[pert.p]*pert.d
                mm.localNF(p,SS)
                SS_nFBp, _ = fn.find_equilibrium(p, SS, mm.nFB)
                nFB_p = mm.out_FB(SS_nFBp)

            # Calculate CoRa
                CoRa = log10(FB_p/FB) / log10(nFB_p/nFB)

        # If not, report error
            else
                println("Unknown error")  
                CoRa = NaN 
                FB = NaN
                error = 3    
            end 
        end
        return CoRa, SS[:, end], FB, error
    end



# --- Generate CoRa curve ---
#     p ----- System parameters
#     u0 ---- Initial conditions
#     mm ---- Models
#     pert -- Perturbation details
#     Returns:
#        CoRa curve,
#        steady states,
#        errors

    function CoRacurve(p, u0, mm, pert)
        r = 10 .^ range(pert.r[1], pert.r[2], length=pert.coras)
        curve = zeros(length(r))
        SSs = zeros(length(r))
        errors = zeros(length(r))
    # Iterate over perturbation range
        for j in 1:length(r)
            p[pert.c] = r[j]
            CoRa = evalCoRa(p, u0, mm, pert)
            curve[j] = CoRa[1]
            u0 =  CoRa[2]
            SSs[j] = CoRa[3]
            errors[j] = CoRa[4]
        end  
        return curve, SSs, errors
    end



# --- Calculate metrics of CoRa curve ---
#     curve -- CoRa curve
#     SSs ---- steady states of controlled variable
#     pert --- Perturbation details
#     Returns:
#        [1] robustness (fraction of points with CoRa<=eps),
#        [2] minRange (minimum rho with CoRa<=eps),
#        [3] maxRange (maximum rho with CoRa<=eps),
#        [4] min(CoRa),
#        [5] optimalRho (rho with min(CoRa)),
#        [6] neg (number of points with negative CoRa),
#        [7] os (number of points with oscillations),
#        [8] e (number of points with other errors),
#        [9] ss (mean steady state of controlled variable for points with CoRa<=eps)

    function metrics(curve, SSs, errors, pert)
        neg = count(x -> x == 1, errors)        
        os = count(x -> x == 2, errors)        #oscillations: number of points with oscillations
        e = count(x -> x == 3, errors)         #other: number of points with other errors

        r = 10 .^ range(pert.r[1], pert.r[2], length=pert.coras)

		i = curve .<= pert.eps;                # indexes less than eps
		j = findall(i);                        # values of rho less than eps
		
        if isempty(j)
			return [0, NaN, NaN, NaN, NaN, neg, os, e, NaN]
        else
            x = copy(curve);
            x[x .=== NaN] .= Inf;                  # replace NaN with Inf for min calculation

            filtered_SSs = SSs[i]                  # steady states for values with CoRa <= eps
            ss = mean(filtered_SSs)                # ss for the values below eps
            return [length(curve[i])/length(curve), r[j[1]], r[j[end]], minimum(x), r[argmin(x)],  neg, os, e, ss]
		end    
    end



# --- Explore ---
#     iARG -- Main arguments
#     mm ---- Models
#     p ----- Core parameters
#     pert -- Perturbation details
#     expl -- Exploration details
#     u0 ---- Initial conditions
#     Returns:
#        Output file with exploration results

    function explore(iARG, mm, p, pert, expl, u0)
        ran = 10 .^ range(pert.r[1], pert.r[2], length=pert.coras)
        n = length(expl.pOp)
    # Open output file
        open(string("./Output/OUT_ExplCoRa_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".txt"), "w") do io
            if (expl.prtD==1)
                    writedlm(io, [vcat([string(param) for param in expl.pOp],"robustness","minRange", "maxRange","min(CoRa)", "optimalRho", "negative_sol", "oscilations", "other_errors", "steady_state", ran)], '\t')
            else
                    writedlm(io, [vcat([string(param) for param in expl.pOp],"robustness","minRange", "maxRange","min(CoRa)", "optimalRho", "negative_sol", "oscilations", "other_errors", "steady_state")], '\t')
            end

        # Set of parameters to explore
            sobol = SobolSeq(n)
            sobol_p = [
                    round.(10.0 .^ (expl.pMin .+ (expl.pMax .- expl.pMin) .* next!(sobol)), digits=4)
                    for _ in 1:expl.n_points
                ]

        # Iterate over set of parameters
            p_orig = copy(p)
            for i in 1:expl.n_points
                println(i)
            # Update parameter values
                for j in 1:n
                    p[expl.pOp[j]] = sobol_p[i][j]
                end
            # Calculate CoRa curve
                curve = CoRacurve(p, u0, mm, pert)

                m = metrics(curve[1], curve[2], curve[3], pert)    # Calculate metrics of curve
                p = copy(p_orig)

            # If printing full CoRa curve specified:
                if (expl.prtD==1)	
                    writedlm(io, [vcat(sobol_p[i], m, curve[1])],'\t')
                else
                    writedlm(io, [vcat(sobol_p[i], m)],'\t')
                end
            end
        end
    end




# --- Optimize ---
#     iARG -- Main arguments
#     mm ---- Models
#     p ----- Core parameters
#     pert -- Perturbation details
#     opt --- Optimization details
#     u0 ---- Initial conditions
#     Returns:
#        Output file with optimization results

    function optimize(iARG, mm, p, pert, opt, u0)
        ran = 10 .^ range(pert.r[1], pert.r[2], length=pert.coras)

        open(string("./Output/OUT_OptCoRa_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".txt"), "w") do io
            if (opt.prtD==1)
                    writedlm(io, [vcat("Iteration", [string(param) for param in opt.pOp],"robustness",string("|CoRa<=",pert.eps,"|"),"min(CoRA)", ran)], '\t')
            else
                    writedlm(io, [vcat("Iteration", [string(param) for param in opt.pOp],"robustness",string("|CoRa<=",pert.eps,"|"),"min(CoRa)")], '\t')
            end

            if opt.rand == 1
                for i in 1:length(opt.pOp)
			    	p[opt.pOp[i]] = round.(10 .^ (rand(Uniform(opt.pMin[i], opt.pMax[i]))), digits = 4);
                end
		    end

            p_copy = copy(p)

            curve0 =  CoRacurve(p, u0, mm, pert)
            m0 = metrics(curve0[1], curve0[2], curve0[3], pert)    # Calculate metrics of curve
            best_rob = m0[1]             # Best robustness
            best_curve = curve0[1]       # Best curve

        # Initial metrics
			op0 = log10(m0[3]/m0[2])    # Property to optimize, range with CoRa<=eps 
			mi0 = m0[4]                 # Minimum CoRa
            r0 = zeros(length(opt.pOp)) 

            if m0[6] != 0 
                println("Negative solutions, try to start with anther parametes")
                return NaN
            end
            if m0[7] != 0 # 
                println("Oscilations, try to start with anther parameters")
                return NaN
            end
            if m0[8] != 0 # there are other errors
                println("Unknown errors, try to start with anther parameters")
                return NaN
            end

            if (opt.prtD==1)	
                writedlm(io, [vcat(0, [p[i] for i in opt.pOp], m0[1], op0, mi0, curve0)],'\t')
            else			
                writedlm(io, [vcat(0, [p[i] for i in opt.pOp], m0[1], op0, mi0)],'\t')
            end

        # Iterate optimization
            p = copy(p_copy)
            for i in 1:opt.iter
                println(i)
            # Generate random changes in parameters according to a normal distribution
                rI = rand(MvNormal(zeros(length(opt.pOp)), zeros(length(opt.pOp)) .+ opt.cov)) 
			# Update parameter values	
                for pI in 1:length(opt.pOp)  
					r0[pI] = p[opt.pOp[pI]]                # Save previous value
					p[opt.pOp[pI]] *= (opt.M .^ rI[pI])    # Update value
                    p[opt.pOp[pI]] = round(p[opt.pOp[pI]]; sigdigits=4)
                    #println("p[",opt.pOp[pI],"] = ", p[opt.pOp[pI]])
				# Exclude values outside regime of exploration:
					if p[opt.pOp[pI]] < (10.0 ^ opt.pMin[pI])
						p[opt.pOp[pI]] = (10.0 ^ opt.pMin[pI])
					elseif p[opt.pOp[pI]] > (10.0 ^ opt.pMax[pI])
						p[opt.pOp[pI]] = (10.0 ^ opt.pMax[pI])
					end
				end

            # Calculate new CoRa curve
                curve1 = CoRacurve(p, u0, mm, pert)
                m1 = metrics(curve1[1], curve1[2], curve1[3], pert)    # Calculate metrics of curve

                op1 = log10(m1[3]/m1[2])  
                mi1 = m1[4];  

            # C1 = minimize min(CoRa) when op0 and op1 = NaN
            # NOTE: As mi0,mi1=[0,1], correct exponential with the expected variance of ~U(0,1)
                xiC = (mi0 ^ 2) / (2 * 0.083)
                xiP = (mi1 ^ 2) / (2 * 0.083)
				c1 = isnan(op0+op1) && (rand() < exp((xiC - xiP)))

            # C2 = maximize range when op0 and op1 != NaN 
			# NOTE: As op0,op1=[0,rrO], but still correct exponential with the expected variance of ~U(0,1)
			#       !! ~U(0,1)*(rrO^2) variance resulted in very noisy runs...
					rrO = pert.r[2] - pert.r[1]
					xiC = (rrO - op0) / (2 * 0.083)
					xiP = (rrO - op1) / (2 * 0.083)
				c2 = rand() < exp((xiC - xiP))

            # C3 = no negative solutions
                c3 = m1[6] == 0
            # C4 = no oscilations
                c4 = m1[7] == 0
            # C5 = no other errors
                c5 = m1[8] == 0 

            # If conditions met to accept new parameter values
				if((c1 || c2 ) && c3 && c4 && c5) 
					op0 = op1
					mi0 = mi1
                    best_rob = m1[1]
                    best_curve = curve1[1]
                    println("ok")
				else
				# If not, revert to previous parameter values
					for pI in 1:length(opt.pOp)
						p[opt.pOp[pI]] = r0[pI]
					end
                    println("no ok")
				end

            # Save results of each iteration
                if (opt.prtD==1)
                    writedlm(io, [vcat(i,[p[x] for x in opt.pOp],best_rob, op0, mi0, best_curve)],'\t')
                else			
                    writedlm(io, [vcat(i,[p[x] for x in opt.pOp],best_rob, op0, mi0)],'\t')
                end
            end
        end
    end



# --- Dynamics ---
#     mm ---- Models
#     u0 ---- Initial conditions
#     p ----- Core parameters
#     tspan - tspan for dynamics
#     pert -- Perturbation details
#     iARG -- Main arguments
#     Returns:
#        Dynamics plot

    function dynamics(mm, u0, p, pert, dyn, iARG)

        p_values = collect(values(p))

        SS,_ = fn.find_equilibrium(p, u0, mm.FB)
        println("Steady state values: ", SS)
        keys_list = collect(keys(p))
        i = findfirst(isequal(pert.p), keys_list)

        prob = ODEProblem(mm.FB, SS, dyn.tspan, p_values)

        dosetimes = dyn.pert_time
        affect!(integrator) = integrator.p[i] *= dyn.pert_size
        cb = PresetTimeCallback(dosetimes, affect!)

        sol = solve(prob, Rodas5Pr(), callback = cb)
        
        plot(sol[dyn.plot,:], xlabel="Time", ylabel="Concentration", label = "Feedback")

        mm.localNF(p,SS)      # Adjust parameters for no feedback system
        p_values = collect(values(p))

        prob = ODEProblem(mm.nFB, SS, dyn.tspan, p_values)
        sol = solve(prob, Rodas5Pr(), callback = cb)
        plot!(sol[dyn.plot,:], linestyle = :dash, label = "No Feedback")
        savefig(string("./Output/OUT_dynamics_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".png"))
    end



    function curve(p, u0, mm, pert, iARG)
        ran = 10 .^ range(pert.r[1], pert.r[2], length=pert.coras)
        c = CoRacurve(p, u0, mm, pert)
        plot(ran, c[1];
            xscale = :log10,
            xlims = (minimum(ran), maximum(ran)),
            ylims = (0, 1),
            xlabel = "P["*string(pert.c)*"]",
            ylabel = "CoRa")
        savefig(string("./Output/OUT_curve_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".png"))
    end

end
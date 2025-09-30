
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


#--- Solve for equilibrium (fast) ---
#    p ------ Set of parameters
#    u0 ----- Initial conditions
#    system - set of ODEs
#    Returns: 
#       vector with equilibrium points

    function find_equilibrium(p, u0, system)
        par = collect(values(p))
        function equilibrium_condition(u)
            du = similar(u)
            system(du, u, par, 0.0)
            return du
        end
        result = nlsolve(equilibrium_condition, u0, autodiff=:forward)
        return result.zero
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

    # Solve Feedback system
        SS = fn.find_equilibrium(p, u0, mm.FB)
        FB = mm.out_FB(SS)  

    # If negative values, try new initial conditions
        i = 1
        while any(x -> x < 0, SS)
            println("negative values in SS, trying new initial conditions")
            u_ = length(u0)
            u_ = fill(10.0^i, u_)

            SS = fn.find_equilibrium(p, u_, mm.FB)
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

    # Calculate CoRa only if no errors
        if error == 0   
        # Solve NF system    
            mm.localNF(p,SS)      # Adjust parameters for no feedback system
            SS_nFB = fn.find_equilibrium(p, SS, mm.nFB)
            nFB = mm.out_FB(SS_nFB)

        # If relative difference between FB and nFB is less than 0.0001, proceed to perturbation   
            if (abs(FB - nFB)/ FB) < 0.0001               
            # perturbation to FB
                p = copy(copy_p)
                p[pert.p] = p[pert.p]*pert.d
                SS_FBp = fn.find_equilibrium(p, SS, mm.FB)
                FB_p = mm.out_FB(SS_FBp)

            # perturbation to nFB    
                p = copy(copy_p)
                p[pert.p] = p[pert.p]*pert.d
                mm.localNF(p,SS)
                SS_nFBp = fn.find_equilibrium(p, SS, mm.nFB)
                nFB_p = mm.out_FB(SS_nFBp)

            # Calculate CoRa
                CoRa = log10(FB_p/FB) / log10(nFB_p/nFB)

        # If not, check for oscillations or other errors
            else
                println("FB and nFB not equal, CoRa not defined")
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
                    println("Unknown error")  
                    CoRa = NaN 
                    FB = NaN
                    error = 3
                end     
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
#        final steady state to use as initial condition in next run

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
        return curve, SSs, u0, errors
    end



# --- Calculate metrics of CoRa curve ---
#    curve -- CoRa curve
#    SSs ---- steady states of controlled variable
#    pert --- Perturbation details
#    Returns:
#       robustness (fraction of points with CoRa<=eps),
#       minRange (minimum rho with CoRa<=eps),
#       maxRange (maximum rho with CoRa<=eps),
#       min(CoRa),
#       optimalRho (rho with min(CoRa)),
#       oscilations (number of points with oscillations),
#       other_errors (number of points with other errors),
#       steady_state (mean steady state of controlled variable for points with CoRa<=eps)

    function metrics(curve, SSs, errors, pert)
        neg = count(x -> x == 1, errors)        
        os = count(x -> x == 2, errors)        #oscillations: number of points with oscillations
        e = count(x -> x == 3, errors)         #other: number of points with other errors

        r = 10 .^ range(pert.r[1], pert.r[2], length=pert.coras)

		i = curve .<= pert.eps;
		j = findall(i);
		
        x = copy(curve);
		x[x .=== NaN] .= Inf;

        indices = curve .< pert.eps          # Condición para filtrar valores en 'a' menores a 0.1
        filtered_SSs = SSs[indices] 
        ss = mean(filtered_SSs)      # ss for the values below 0.1
    
		if isempty(j)
			return [length(curve[i])/length(curve), NaN, NaN, minimum(x), r[argmin(x)], os, other, NaN]
		end

		return [length(curve[i])/length(curve), r[j[1]], r[j[end]], minimum(x), r[argmin(x)],  os, other, ss]
    end




    function explore(iARG, mm, p, pert, expl, u0)
        ran = 10 .^ range(pert.r[1], pert.r[2], length=pert.coras)

        open(string("./output/OUT_ExplCoRa_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".txt"), "w") do io
            if (expl.prtD==1)
                    writedlm(io, [vcat([string(param) for param in expl.pOp],"robustness","minRange", "maxRange","min(CoRa)", "optimalRho", "oscilations", "other_errors", "steady_state", ran)], '\t')
            else
                    writedlm(io, [vcat([string(param) for param in expl.pOp],"robustness","minRange", "maxRange","min(CoRa)", "optimalRho", "oscilations", "other_errors", "steady_state")], '\t')
            end

            # Set of parameters    
            sobol = SobolSeq(expl.n_params)
            sobol_p = [
                    round.(10.0 .^ (expl.pMin .+ (expl.pMax .- expl.pMin) .* next!(sobol)), digits=4)
                    for _ in 1:expl.n_points
                ]

            p_orig = copy(p)
            for i in 1:expl.n_points
                println(i)
                for j in 1:expl.n_params
                    p[expl.pOp[j]] = sobol_p[i][j];
                end
                curve, SSs, u1 = CoRacurve(p, u0, mm, pert)

                #println(u1)
                m = metrics(curve, SSs, pert);    # Calculate metrics of DY curve
               
                p = copy(p_orig)

                if (expl.prtD==1)	# If printing full CoRa curve specified:
                    writedlm(io, [vcat(sobol_p[i], m, curve)],'\t')
                else			# Else:
                    writedlm(io, [vcat(sobol_p[i], m)],'\t')
                end

            end
        end
    end



    function optimize(iARG, mm, p, pert, opt, u0)
        ran = 10 .^ range(pert.r[1], pert.r[2], length=pert.coras)

        open(string("./output/OUT_OptCoRa_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".txt"), "w") do io
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

            curve0, SSs0, u1 = CoRacurve(p, u0, mm, pert)
            m0 = metrics(curve0, SSs0, pert);    # Calculate metrics of DY curve

			op0 = log10(m0[3]/m0[2]);   # Property to optimize (e.g. DY<=eps range length)
			mi0 = m0[4];  
            r0 = zeros(length(opt.pOp));

            if m0[6] != 0 # there are oscillations
                println("oscilations, try to start with anther parametes")
                return NaN
            end
            if m0[7] != 0 # there are other errors
                println("there are other errors, try to start with anther parameters")
                return NaN
            end

            if (opt.prtD==1)	
                writedlm(io, [vcat(0, [p[i] for i in opt.pOp], m0[1],op0, mi0, curve0)],'\t')
            else			
                writedlm(io, [vcat(0, [p[i] for i in opt.pOp], m0[1], op0, mi0)],'\t')
            end


            p = copy(p_copy)
            for i in 1:opt.iter
                println(i)
                rI = rand(MvNormal(zeros(length(opt.pOp)), zeros(length(opt.pOp)) .+ opt.cov)); # Random values to update parameters
				for pI in 1:length(opt.pOp)  # Update parameter values
					r0[pI] = p[opt.pOp[pI]]; # Save previous value
					p[opt.pOp[pI]] *= (opt.M .^ rI[pI]); # Update value
                    p[opt.pOp[pI]] = round(p[opt.pOp[pI]]; sigdigits=4)
                    println("p[",opt.pOp[pI],"] = ", p[opt.pOp[pI]])
					# Exclude values outside regime of exploration:
					if p[opt.pOp[pI]] < (10.0 ^ opt.pMin[pI])
						p[opt.pOp[pI]] = (10.0 ^ opt.pMin[pI])
					elseif p[opt.pOp[pI]] > (10.0 ^ opt.pMax[pI])
						p[opt.pOp[pI]] = (10.0 ^ opt.pMax[pI])
					end
				end
                
                curve1, SSs1, u1 = CoRacurve(p, u0, mm, pert)
                m1 = metrics(curve1, SSs1, pert);    # Calculate metrics of DY curve

                op1 = log10(m1[3]/m1[2]);   # Property to optimize (e.g. DY<=eps range length)
                mi1 = m1[4];  

                # Evaluate if accept new parameter values or not:
				## Only accept in the regime of interest, i.e. DY>=0:
                c1 = (mi1>=0);
				## If DY>eps for all rho, evaluate the min(DY) for both sets:
				### NOTE: As mi0,mi1=[0,1], correct exponential with the expected variance of ~U(0,1)
					xiC = (mi0 ^ 2) / (2 * 0.083);
					xiP = (mi1 ^ 2) / (2 * 0.083);
				c2 = isnan(op0+op1) && (rand() < exp((xiC - xiP)));
				## If DY>=eps for some rho, evaluate the |DY<=eps| for both sets:
				### NOTE: As op0,op1=[0,rrO], but still correct exponential with the expected variance of ~U(0,1)
				###       !! ~U(0,1)*(rrO^2) variance resulted in very noisy runs...
					rrO = pert.r[2] - pert.r[1];
					xiC = (rrO - op0) / (2 * 0.083);
					xiP = (rrO - op1) / (2 * 0.083);
				c3 = rand() < exp((xiC - xiP));

                c4 = m1[6] == 0 # there are not oscillations
                c5 = m1[7] == 0 # there are not other errors

				if(c1 && (c2 || c3) && c4 && c5) 
					# If yes, update "reference" system
					op0 = op1;
					mi0 = mi1;
                    m0[1] = m1[1];
					curve0 = curve1;
                    println("ok")
				else
					# If not, revert to previous parameter values
					for pI in 1:length(opt.pOp)
						p[opt.pOp[pI]] = r0[pI];
					end
                    println("no ok")
				end

               
                if (opt.prtD==1)	# If printing full CoRa curve specified:
                    writedlm(io, [vcat(i,[p[x] for x in opt.pOp],m0[1], op0, mi0, curve0)],'\t')
                else			
                    writedlm(io, [vcat(i,[p[x] for x in opt.pOp],m0[1], op0, mi0)],'\t')
                end
            end
        end
    end

    function dynamics(mm, u0, p, tspan, pert, iARG)
        p[:dyn] = 0
        p_values = collect(values(p))
        tspan = (0.0, 500000.0)
        #println(evalCoRa(p, u0, mm, pert))

        prob = ODEProblem(mm.FB, u0, tspan, p_values)

        sol = solve(prob, Rodas5())

        SS = sol[1,end]
        println("Generating dynamics plot...", SS)

		p[:miss] = sol[4,end]  


        p[:dyn] = 1
        p_values = collect(values(p))
        tspan = (0.0, 10000.0)
        prob = ODEProblem(mm.FB, sol[:,end], tspan, p_values)
        sol = solve(prob, Rodas5())  # Usar un solver robusto

        t = sol.t
        y = sol[1, :]   # primera variable, cambia el índice si quieres otra
        # filtrar solo los últimos tiempos
        mask = t .>= 0
        plot(t[mask], y[mask], xlabel="Time", ylabel="Concentration")

        prob = ODEProblem(mm.nFB, sol[:,end], tspan, p_values)
        sol = solve(prob, Rodas5())  # Usar un solver robusto
        t = sol.t
        y = sol[1, :]   # primera variable, cambia el índice si quieres otra
        # filtrar solo los últimos tiempos
        mask = t .>= 0
        #plot!(t[mask], y[mask], xlabel="Time", ylabel="Concentration")

        savefig(string("OUT_Dyn_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".png"))
        
        println(evalCoRa(p, u0, mm, pert))

        p[:dyn] = 0
    end

    function curve(p, u0, mm, pert, iARG)
        ran = 10 .^ range(pert.r[1], pert.r[2], length=pert.coras)
        curve, SSs, u1 = CoRacurve(p, u0, mm, pert)
        plot(ran, curve;
            xscale = :log10,
            xlims = (minimum(ran), maximum(ran)),
            ylims = (0, 1),
            xlabel = "P",
            ylabel = "CoRa")
        savefig(string("OUT_curve_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".png"))
    end

end


    #=
    #--- Solve for equilibrium (slow) ---
    # p ------ Set of parameters
    # u0 ----- Initial conditions
    # system - set of ODEs
    # Returns: vector with equilibrium points 

    function solve_to_steady_state(p, u0, system, u0, p, tspan)
        p_values = collect(values(p))
        prob = ODEProblem(system!, u0, tspan, p_values)
        sol = solve(prob, Rodas5(), reltol=1e-8, abstol=1e-8)  # Usar un solver robusto
        # Verificar si el sistema alcanzó el estado estacionario
        du = similar(u0)
        system!(du, sol.u[end], p_values, sol.t[end])
        if maximum(abs.(du)) < 1e-6
            #println("El sistema alcanzó el estado estacionario.")
            return sol.u[end]
        else
            println("no ss")
            return 0.0
        end
    end




        function Check(ssR, soR, rtol)
        if isnan(ssR) || isnan(soR) || (abs(ssR) - abs(soR)) > 1e-4
            rtol *= 1e-3
            if(rtol < 1e-24)
                println("ERROR: Check NF system (reltol=",rtol,").")
                #println(vcat(pert.p,i,[p[eval(Meta.parse(string(":",i)))] for i in syst.sys.ps],mm.outFB(ssR),mm.outNF(soR)))
                if(abs(ssR - soR)/ssR > 0.01)
                    ssR, soR = Restart(ssR, soR)
                    println("Error too large. SS results excluded!")
                end
            end
            return ssR, soR, rtol, "Insufficient"
        else
            return ssR, soR, rtol, "Sufficient"
        end
    end;



        
    function solver(type, p, u0, mm,  model, pert)
        if type == "fast"
            SS = fn.find_equilibrium(p, u0, model)
            FB = mm.outFB_fast(SS)
        else
            SS = fn.solve_to_steady_state(model, u0, p, pert.tspan)
            FB = mm.outFB_slow(SS)
        end
        return SS, FB
    end
    =#
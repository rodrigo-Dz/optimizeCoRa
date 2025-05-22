module fn
	# Required libraries
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



        # Solve 
    function solve_to_steady_state(system!, u0, p, tspan)
        p_values = collect(values(p))
        prob = ODEProblem(system!, u0, tspan, p_values)

        sol = solve(prob, Rodas5(), reltol=1e-8, abstol=1e-8)  # Usar un solver robusto

        # Verificar si el sistema alcanzó el estado estacionario
        du = similar(u0)
        system!(du, sol.u[end], p_values, sol.t[end])
        if maximum(abs.(du)) < 1e-6
            #println("El sistema alcanzó el estado estacionario.")
            return sol
        else
            println(" no ss")
            return 0.0
        end

        
    end

    # Encontrar puntos de equilibrio para el primer sistema
    function find_equilibrium(p, u0, system!)
        p_values = collect(values(p))
        function equilibrium_condition(u)
            du = similar(u)
            system!(du, u, p_values, 0.0)
            return du
        end
        result = nlsolve(equilibrium_condition, u0, autodiff=:forward)
        return result.zero
    end


    # Calcular el Jacobiano en un punto de equilibrio
    function compute_jacobian(u_eq, p, system!)
        p_values = collect(values(p))
        function equilibrium_condition(u)
            du = similar(u)
            system!(du, u, p_values, 0.0)
            return du
        end
        ForwardDiff.jacobian(equilibrium_condition, u_eq)
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


    function evalCoRa(p, u0, mm, pert)
        copy_p = copy(p)
        # SS FB system
        if pert.solver == "fast"
            SS = fn.find_equilibrium(p, u0, mm.FB)
            FB = mm.outFB_fast(SS)
        else
            SS = fn.solve_to_steady_state(p, u0, mm.FB, pert.tspan)
            FB = mm.outFB_slow(SS)
        end

        p = copy(copy_p)
        mm.localNF(p,SS)

        # SS NF system
        if pert.solver == "fast"
            SS_nFB = fn.find_equilibrium(p, SS, mm.nFB)
            nFB = mm.out_nFB_fast(SS_nFB)
        else
            SS_nFB = fn.solve_to_steady_state(p, SS, mm.nFB, pert.tspan)
            nFB = mm.out_nFB_slow(SS_nFB)
        end

        if (abs(FB - nFB)/ max(abs(FB),abs(nFB))) < 0.0001
            p = copy(copy_p)
            p[pert.p] = p[pert.p]*pert.d
            
            if pert.solver == "fast"
                SS_FBp = fn.find_equilibrium(p, SS, mm.FB)
                FB_p = mm.out_nFB_fast(SS_FBp)
            else
                SS_FBp = solve_to_steady_state(p, SS, mm.FB, pert.tspan)
                FB_p = mm.out_nFB_slow(SS_FBp)
            end

            p = copy(copy_p)
            p[pert.p] = p[pert.p]*pert.d
            mm.localNF(p,SS)


            if pert.solver == "fast"
                SS_nFBp = fn.find_equilibrium(p, SS, mm.nFB)
                nFB_p = mm.out_nFB_fast(SS_nFBp)
            else
                SS_nFBp = fn.solve_to_steady_state(p, SS, mm.nFB, pert.tspan)
                nFB_p = mm.out_nFB_slow(SS_nFBp)
            end

            if (FB_p < 0) || (FB<0) || (nFB_p<0) || (nFB<0)
                println("error en uno de estos:", FB, FB_p, nFB_p, nFB)
                CoRa  = NaN
                SS_controled = NaN
            else
                CoRa = log10(FB_p/FB) / log10(nFB_p/nFB)
                SS_controled = FB
            end
        
        else
            J = fn.compute_jacobian(SS, p, mm.FB)
            eigenvalues = eigvals(J)  # Calcular autovalores
            has_oscilations = any(imag(l) != 0 && real(l) > 0 for l in eigenvalues) # Eigenvalues complejos con parte real positiva
            if has_oscilations == true
                CoRa = 2  
                SS_controled = NaN 
                println("oscila")
            else
                println("No oscila pero pasa algo raro, intenta con otro solver")  # Todo lo demás
                CoRa = 3  #mark other type of erros
                SS_controled = NaN
            end     

        end 

        return CoRa, SS, SS_controled
    end

    function CoRacurve(p, u0, mm, pert)
        ran = 10 .^ range(pert.r[1], pert.r[2], length=pert.coras)
        curve = zeros(length(ran))
        SSs = zeros(length(ran))
        for j in 1:length(ran)
            p[pert.c] = ran[j]
            CoRa_r = evalCoRa(p, u0, mm, pert)
            curve[j] = CoRa_r[1]
            SSs[j] = CoRa_r[3]
            u0 =  CoRa_r[2]
        end  
        return curve, SSs
    end


    function metrics(curve, SSs, pert)
        r = 10 .^ range(pert.r[1], pert.r[2], length=pert.coras)
		i = curve .<= pert.eps;
		j = findall(i);
		x = copy(curve);
		x[x .=== NaN] .= Inf;

        os = count(x -> x == 2, curve)        #oscillations: number of points with oscillations
        other = count(x -> x == 3, curve)    #other: number of points with other errors
        indices = curve .< pert.eps          # Condición para filtrar valores en 'a' menores a 0.1
        filtered_SSs = SSs[indices] 
        ss = mean(filtered_SSs)      # ss for the values below 0.1
    
		if isempty(j)
			return [length(curve[i])/length(curve), NaN, NaN, minimum(x), r[argmin(x)], os, other, NaN]
		end

		return [length(curve[i])/length(curve), r[j[1]], r[j[end]], minimum(x), r[argmin(x)],  os, other, ss]
    end



    function explore_all(iARG, mm, p, pert, expl, u0)
        ran = 10 .^ range(pert.r[1], pert.r[2], length=pert.coras)

        open(string("./output/OUT_ExplCoRa_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".txt"), "w") do io
            if (expl.prtD==1)
                    writedlm(io, [vcat([string(param) for param in expl.pOp],"proportion<=eps","minRange", "maxRange","min(CoRa)", "optimalRho", "oscilations", "other_errors", "steady_state", ran)], '\t')
            else
                    writedlm(io, [vcat([string(param) for param in expl.pOp],"proportion<=eps","minRange", "maxRange","min(CoRa)", "optimalRho", "oscilations", "other_errors", "steady_state")], '\t')
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

                curve, SSs = CoRacurve(p, u0, mm, pert)

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
                    writedlm(io, [vcat("Iteration", [string(param) for param in opt.pOp],"proportion<=eps",string("|CoRa<=",pert.eps,"|"),"min(CoRA)", ran)], '\t')
            else
                    writedlm(io, [vcat("Iteration", [string(param) for param in opt.pOp],"proportion<=eps",string("|CoRa<=",pert.eps,"|"),"min(CoRa)")], '\t')
            end

            if opt.rand == 1
                for i in 1:length(opt.pOp)
			    	p[opt.pOp[i]] = round.(10 .^ (rand(Uniform(opt.pMin[i], opt.pMax[i]))), digits = 4);
                end
		    end

            p_copy = copy(p)

            curve0, SSs0 = CoRacurve(p, u0, mm, pert)
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
					# Exclude values outside regime of exploration:
					if p[opt.pOp[pI]] < (10.0 ^ opt.pMin[pI])
						p[opt.pOp[pI]] = (10.0 ^ opt.pMin[pI])
					elseif p[opt.pOp[pI]] > (10.0 ^ opt.pMax[pI])
						p[opt.pOp[pI]] = (10.0 ^ opt.pMax[pI])
					end
				end
                
                curve1, SSs1 = CoRacurve(p, u0, mm, pert)
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


end

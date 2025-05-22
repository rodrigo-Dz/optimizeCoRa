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


    function explore_all(iARG, mm, p, pert, opt, c1)
        ran = 10 .^ range(pert.r[1], pert.r[2], length=pert.coras)

        open(string("./output/OUT_OptCoRa_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".txt"), "w") do io
            if (opt.prtD==1)
                    writedlm(io, [vcat([string(param) for param in opt.pOp],string("|CoRa<=",pert.eps,"|"),"min(CoRA)", "steady_state","oscilations", "other_errors", ran)], '\t')
            else
                    writedlm(io, [vcat([string(param) for param in opt.pOp],string("|CoRa<=",pert.eps,"|"),"min(CoRa)", "steady_state","oscilations", "other_errors")], '\t')
            end

            # Set of parameters    
            sobol = SobolSeq(opt.n_params)
            sobol_p = [
                    round.(10.0 .^ (opt.pMin .+ (opt.pMax .- opt.pMin) .* next!(sobol)), digits=4)
                    for _ in 1:opt.n_points
                ]

            p_orig = copy(p)
            for i in 1:opt.n_points
                for j in 1:opt.n_params
                    p[opt.pOp[j]] = sobol_p[i][j];
                end

                curve = zeros(length(ran))
                SSs = zeros(length(ran))
                u0 = c1

                for j in 1:length(ran)
                    p[pert.c] = ran[j]
                    CoRa_r = evalCoRa(p, u0, mm, pert)

                    curve[j] = CoRa_r[1]
                    SSs[j] = CoRa_r[3]
                    u0 =  CoRa_r[2]
                end  

                p = copy(p_orig)

                println(curve)
                rob = count(x -> x < pert.eps, curve) #robustness: number of points below 0.1
                os = count(x -> x == 2, curve)        #oscillations: number of points with oscillations
                other = count(x -> x == 3, curve)   #other: number of points with other errors

                min = minimum(curve) #min CoRa
                indices = curve .< pert.eps  # Condición para filtrar valores en 'a' menores a 0.1
                filtered_SSs = SSs[indices]
                ss = mean(filtered_SSs)      # ss for the values below 0.1
                
                if (opt.prtD==1)	# If printing full CoRa curve specified:
                    writedlm(io, [vcat(sobol_p[i], rob, min ,ss, os, other, curve)],'\t')
                else			# Else:
                    writedlm(io, [vcat(sobol_p[i], rob, min ,ss, os, other, curve)],'\t')
                end
            end
        end
    end



    function optimize(iARG, mm, p, pert, opt, u0)
        ran = 10 .^ range(pert.r[1], pert.r[2], length=pert.coras)

        open(string("./output/OUT_OptCoRa_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".txt"), "w") do io
            if (opt.prtD==1)
                    writedlm(io, [vcat([string(param) for param in opt.pOp],string("|CoRa<=",pert.eps,"|"),"min(CoRA)", "steady_state","oscilations", "other_errors", ran)], '\t')
            else
                    writedlm(io, [vcat([string(param) for param in opt.pOp],string("|CoRa<=",pert.eps,"|"),"min(CoRa)", "steady_state","oscilations", "other_errors")], '\t')
            end
            #n_points, a, b, n_params, u0,  tspan,  solver, eps)
            results = []
            p_ori = copy(p)

            # Set of parameters    
            sobol = SobolSeq(opt.n_params)
            sobol_p = [
                    10.0 .^ (opt.pMin .+ (opt.pMax .- opt.pMin) .* next!(sobol))
                    for _ in 1:opt.n_points
                ]

            for i in 1:opt.n_points
                p = copy(p_ori)
                for j in 1:opt.n_params
                    p[opt.pOp[j]] = sobol_p[i][j];
                end

                curve = zeros(length(ran))
                SSs = zeros(length(ran))
                u0 = u0
                for j in 1:length(ran)
                    p[pert.c] = ran[j]
                    copy_p = copy(p)  

                    # SS FB system
                    if opt.solver == "fast"
                        SS = fn.find_equilibrium(p, u0, mm.FB)
                        FB = mm.outFB_fast(SS)
                        println(FB)
                    else
                        SS = fn.solve_to_steady_state(p, u0, mm.FB, opt.tspan)
                        FB = mm.outFB_slow(SS)
                    end

                    u0 = SS

                    p = copy(copy_p)
                    mm.localNF(p,SS)

                    if opt.solver == "fast"
                        SS_nFB = fn.find_equilibrium(p, u0, mm.nFB)
                        nFB = mm.out_nFB_fast(SS_nFB)
                    else
                        SS_nFB = fn.solve_to_steady_state(p, u0, mm.nFB, opt.tspan)
                        nFB = mm.out_nFB_slow(SS_nFB)
                    end

                    if (abs(FB - nFB)/ max(abs(FB),abs(nFB))) < 0.0001
                        p = copy(copy_p)
                        p[pert.p] = p[pert.p]*pert.d
                        
                        if opt.solver == "fast"
                            SS_FBp = fn.find_equilibrium(p, u0, mm.FB)
                            FB_p = mm.out_nFB_fast(SS_FBp)
                        else
                            SS_FBp = solve_to_steady_state(p, u0, mm.FB, opt.tspan)
                            FB_p = mm.out_nFB_slow(SS_FBp)
                        end

                        p = copy(copy_p)
                        p[pert.p] = p[pert.p]*pert.d
                        mm.localNF(p,SS)


                        if opt.solver == "fast"
                            SS_nFBp = fn.find_equilibrium(p, u0, mm.nFB)
                            nFB_p = mm.out_nFB_fast(SS_nFBp)
                        else
                            SS_nFBp = fn.solve_to_steady_state(p, u0, mm.nFB, opt.tspan)
                            nFB_p = mm.out_nFB_slow(SS_nFBp)
                        end

                        if (FB_p < 0) || (FB<0) || (nFB_p<0) || (nFB<0)
                            println("error en uno de estos:", FB, FB_p, nFB_p, nFB)
                            curve[j] = NaN
                            SSs[j] = NaN
                        else
                            curve[j] = log10(FB_p/FB) / log10(nFB_p/nFB)
                            SSs[j] = FB
                        end
                    
                    else
                        J = fn.compute_jacobian(SS, p, mm.FB)
                        eigenvalues = eigvals(J)  # Calcular autovalores
                        has_oscilations = any(imag(l) != 0 && real(l) > 0 for l in eigenvalues) # Eigenvalues complejos con parte real positiva
                        if has_oscilations == true
                            curve[j] = 2  
                            SSs[j] = NaN 
                            println("oscila")
                        else
                            println("No oscila pero pasa algo raro, intenta con otro solver")  # Todo lo demás
                            curve[j] = 3  #mark other type of erros
                            SSs[j] = NaN
                        end     

                    end            
                end  
                rob = count(x -> x < pert.eps, curve) #robustness: number of points below 0.1
                os = count(x -> x == 2, curve)        #oscillations: number of points with oscillations
                other = count(x -> x == 3, curve)   #other: number of points with other errors

                min = minimum(curve) #min CoRa
                indices = curve .< pert.eps  # Condición para filtrar valores en 'a' menores a 0.1
                filtered_SSs = SSs[indices]
                ss = mean(filtered_SSs)      # ss for the values below 0.1
                
                if (opt.prtD==1)	# If printing full CoRa curve specified:
                    writedlm(io, [vcat(sobol_p[i], rob, min ,ss, os, other, curve)],'\t')
                else			# Else:
                    writedlm(io, [vcat(sobol_p[i], rob, min ,ss, os, other, curve)],'\t')
                end
            end
        end
    end

end

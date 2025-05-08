## CoRa optimization

using DifferentialEquations
using ForwardDiff
using LinearAlgebra
using Plots
using NLsolve
using Sobol
using DelimitedFiles;
using Statistics
using DataStructures



function explore_all(n_points, a, b, n_params, u0, system_FB, system_nFB, tspan, pOp, p, mY_range, solver)
    results = []
    p_ori = copy(p)
    # Set of parameters    
    sobol_p = SobolSeq(n_params)
    sobol_p = [10.0 .^ (a .+ (b - a) .* next!(sobol_p)) for _ in 1:n_points]

    for i in 1:n_points
        println(i)
        p = copy(p_ori)
        for j in 1:length(pOp)
			#p[mrw.pOp[i]] = 10.0 .^ (rand(Uniform(mrw.pMin[i], mrw.pMax[i])));
            p[pOp[j]] = sobol_p[i][j];

		end

        curve = zeros(length(mY_range))
        SSs = zeros(length(mY_range))
        next_ci = u0
        for j in 1:length(mY_range)
            p[:mY] = mY_range[j]
            copy_p = copy(p)  
            #(sol)
            if solver == "fast"
                SS = find_equilibrium(p, next_ci, system_FB)
                FB = outFB_fast(SS)
            else
                SS = solve_to_steady_state(p, next_ci, system_FB, tspan)
                FB = outFB_slow(SS)
                SS = sol.u[end]
            end
            #sol = solve_to_steady_state(system_FB, next_ci, p, tspan)
            #SS = sol.u[end]  # Usar el último punto de la solución como estado estacionario
            
        
            FB = (SS[1])  # Almacenar el valor de Y en el estado estacionario
            next_ci = SS

            p = copy(copy_p)
            p[:mUs] = p[:mU] * FB

            SS_ = find_equilibrium(p, SS, system_nFB)
            nFB = (SS_[1])


            #sol = solve_to_steady_state(system_nFB, SS, p, tspan)
            #ss_ = sol.u[end]

            if abs(FB - nFB) < 0.0001

                p = copy(copy_p)
                p[:mY] = p[:mY] * 1.05
                sol_ = find_equilibrium(p, SS, system_FB)
                ss_= sol_[1]
                #sol = solve_to_steady_state(system_FB, SS, p, tspan)
                #ss_ = sol.u[end]
                FB_p = (ss_[1])

                p = copy(copy_p)
                p[:mY] = p[:mY] * 1.05
                p[:mUs] = p[:mU] * FB
                sol_ = find_equilibrium(p, SS, system_nFB)
                ss_= sol_[1]
                #sol = solve_to_steady_state(system_nFB, SS, p, tspan)
                #ss_ = sol.u[end]
                nFB_p = (ss_[1])

                if (FB_p < 0) || (FB<0) || (nFB_p<0) || (nFB<0)
                    println("error en uno de estos:", FB, FB_p, nFB_p, nFB)
                    curve[j] = NaN
                    SSs[j] = NaN
                else
                    curve[j] = log10(FB_p/FB) / log10(nFB_p/nFB)
                    SSs[j] = FB
                end
            else
                J = compute_jacobian(sol, p, system_FB)
                eigenvalues = eigvals(J)  # Calcular autovalores
                has_oscilations = any(imag(l) != 0 && real(l) > 0 for l in eigenvalues)
                if has_oscilations == true
                    curve[j] = 1
                    SSs[j] = NaN # Eigenvalues complejos con parte real positiva
                    println("oscila")
                else
                    println("No oscila pero pasa algo raro, intenta con otro solver")  # Todo lo demás
                    curve[j] = 0
                    SSs[j] = NaN
                end     
                println("no son iguales")
                #curve[j] = NaN
                #SSs[j] = NaN
            end            
        end  
        r_01 = count(x -> x < 0.1, curve)
        r_02 = count(x -> x < 0.2, curve)
        r_03 = count(x -> x < 0.3, curve)
        os = count(x -> x == 1, curve)

        indices = curve .< 0.1  # Condición para filtrar valores en 'a' menores a 0.1
        # Obtener los valores filtrados de 'b'
        filtered_SSs = SSs[indices]
        # Calcular el promedio de los valores filtrados de 'b'
        ss = mean(filtered_SSs)
        push!(results,  vcat(collect(values(p)), r_01, r_02, r_03, os, ss))
    end

    # Convertir results a una matriz
    results_matrix = reduce(vcat, transpose.(results))
    # Guardar los resultados en un archivo CSV
    writedlm("bifurcation_results_ATFv1_f.csv", results_matrix, ',')

    return results_matrix

end



using Pkg

Pkg.activate("../.")

using DifferentialEquations
using Random, Statistics, Distributions
using Optim, DelimitedFiles

DiffEq = DifferentialEquations

include("temp_assembly.jl")

Random.seed!(1)
N = 100
x0 = rand(N)
tspan = (0, 1e6)

#simulation parameters
N_T, N_σ, N_rep = 25,3,5
T_vec = range(280,300,length = N_T)
σ_vec = [0.01, 0.1, 0.25]

#Temperature dependencies
function r_func_temp(T, E_σ, x...)
    return(boltz.(rand(LogNormal(0.0,0.5), x), rand(Normal(0.2, E_σ), x), KtoT(T, 290.0)))
end

function int_func_temp(T, E_σ, p_neg, x...)
    a = boltz.(rand(LogNormal(-0.1,0.05), x), rand(Normal(0.2,E_σ), x), KtoT(T, 290.0))
#     a[rand(Bool, x)] .*= -1.0
    return(-a ./ 100)
end

#results array
Nsp_res = Array{Any,3}(undef, N_T, N_rep, N_σ)
t_res = Array{Any,3}(undef, N_T, N_rep, N_σ)
params_mat = Array{Any,3}(undef,N_T,N_rep, N_σ)

for i = 1:N_T
    Threads.@threads for j = 1:N_rep
        for k = 1:N_σ
            @show (T_vec[i] , j, k, Threads.threadid())
            #define functions 
            r_func(x...) = r_func_temp(T_vec[i], σ_vec[k], x...)
            #0.6 of interactions are competitive
            int_func(x...) = int_func_temp(T_vec[i], σ_vec[k], 0.6, x...)

            #generate params
            r = r_func(N)
            a = int_func(N,N)

            #set intra-specific interactions
            [a[x,x] = -1.0 for x = 1:N]


            p = Param(N, r, a, 0.0)

            prob = DiffEq.ODEProblem(dx!, x0, tspan, p)
            sol = DiffEq.solve(prob, save_everystep = false)
            @show( sum(sol[end] .> eps()) )
            
            #save solution
            Nsp_res[i,j,k] = string(sum(sol[end] .> eps()))
            t_res[i,j,k] = string(sol.t)
            params_mat[i,j,k] = string( (mean(int_func(10000)), std(r_func(10000))) ) 
        end
    end
end

writedlm("HelperScripts/Data/Non_assembly/N_1.csv", Nsp_res[:,:,1], ",")
writedlm("HelperScripts/Data/Non_assembly/N_2.csv", Nsp_res[:,:,2], ",")
writedlm("HelperScripts/Data/Non_assembly/N_3.csv", Nsp_res[:,:,3], ",")

writedlm("HelperScripts/Data/Non_assembly/t_1.csv", t_res[:,:,1], ",")
writedlm("HelperScripts/Data/Non_assembly/t_2.csv", t_res[:,:,2], ",")
writedlm("HelperScripts/Data/Non_assembly/t_3.csv", t_res[:,:,3], ",")

writedlm("HelperScripts/Data/Non_assembly/params_1.csv", params_mat[:,:,1], ",")
writedlm("HelperScripts/Data/Non_assembly/params_2.csv", params_mat[:,:,2], ",")
writedlm("HelperScripts/Data/Non_assembly/params_3.csv", params_mat[:,:,3], ",")
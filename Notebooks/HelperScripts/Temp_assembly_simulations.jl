using Pkg

Pkg.activate("../.")

using Plots, DifferentialEquations
using Random, Statistics, Distributions, StatsPlots
using Optim, DelimitedFiles

Plots.gr()

DiffEq = DifferentialEquations

include("../HelperScripts/temp_assembly.jl")


Random.seed!(1)
N = 3
x0 = rand(N)
tspan = (0, 1e6)

#simulation parameters
N_T, N_σ, N_rep = 25,3,5
T_vec = range(280,300,length = N_T)
σ_vec = [0.01, 0.05, 0.1]

#Temperature dependencies
r_func_temp(T, E_σ, x...) = boltz.(rand(LogNormal(0.0,0.5), x), rand(Normal(0.0, 0.28), x), KtoT(T, 300.0))
int_func_temp(T, E_σ, x...) = boltz.(rand(LogNormal(-1,0.05), x), rand(Normal(0.65,E_σ), x), KtoT(T, 300.0))

#callback

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
        int_func(x...) = int_func_temp(T_vec[i], σ_vec[k], x...)

        cb = add_at_equi(int_func, r_func, abstol = 1e-6, reltol = 1e-6)

        #generate params
        r = r_func(N)
        a = int_func(N,N)
        #set sign of interactions with p = 0.6
        a[rand(N,N) .<= 0.6] .*= -1    
        #set intra-specific interactions
        [a[x,x] = -1.0 for x = 1:N]
            
        show(a)
        show(mean(r_func(1000)))
        show(std(r_func(1000)))

        p = Param(N, r, a, 0.0)

        prob = DiffEq.ODEProblem(dx!, x0, tspan, p)
        sol = DiffEq.solve(prob, callback = cb, save_everystep = false)

        #save solution
        Nsp_res[i,j,k] = string(get_N(sol))
        t_res[i,j,k] = string(sol.t)
        params_mat[i,j,k] = string( (mean(int_func(10000)), std(r_func(10000))) ) 
        end
    end
end


writedlm("HelperScripts/Data/N_1.csv", Nsp_res[:,:,1], ",")
writedlm("HelperScripts/Data/N_2.csv", Nsp_res[:,:,2], ",")
writedlm("HelperScripts/Data/N_3.csv", Nsp_res[:,:,3], ",")

writedlm("HelperScripts/Data/t_1.csv", t_res[:,:,1], ",")
writedlm("HelperScripts/Data/t_2.csv", t_res[:,:,2], ",")
writedlm("HelperScripts/Data/t_3.csv", t_res[:,:,3], ",")

writedlm("HelperScripts/Data/params_1.csv", params_mat[:,:,1], ",")
writedlm("HelperScripts/Data/params_2.csv", params_mat[:,:,2], ",")
writedlm("HelperScripts/Data/params_3.csv", params_mat[:,:,3], ",")
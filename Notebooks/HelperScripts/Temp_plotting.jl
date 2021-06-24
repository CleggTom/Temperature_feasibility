using Pkg

Pkg.activate(".")

using Plots
using DelimitedFiles
Plots.gr()

#read files
#simulations are done in "./HelperScripts/Temp_assembly_simulations.jl"

N_mat = Array{Any,3}(undef, 25, 5, 3)
t_mat = Array{Any,3}(undef, 25, 5, 3)
p_mat = Array{Any,3}(undef, 25, 5, 3)

for i = 1:3
    #parse N_sp
    N_mat_temp = readdlm( join( ["./Notebooks/HelperScripts/Data/N_",i,".csv"] ), ',' )
    N_mat_temp = Meta.parse.(N_mat_temp) .|> eval
    N_mat[:,:,i] .= N_mat_temp
    
    #parse t
    t_mat_temp = readdlm( join( ["./Notebooks/HelperScripts/Data/N_",i,".csv"] ), ',' )
    t_mat_temp = Meta.parse.(t_mat_temp) .|> eval
    t_mat[:,:,i] .= t_mat_temp
    
    #parse r
    params_mat_temp =  readdlm( join( ["./Notebooks/HelperScripts/Data/N_",i,".csv"] ), ',' )
    params_mat_temp = Meta.parse.(params_mat_temp) .|> eval
    p_mat[:,:,i] .= params_mat_temp
end

p = plot(legend = false)

for i = eachindex(t_mat[:,:,1])
    plot!(t_mat[i], N_mat[i])
end

Plots.savefig(p, "./data/N_t.png")


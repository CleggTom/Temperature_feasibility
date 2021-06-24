####
# Temperature assembly simulations
####

#temperature conversion
KtoT(K, Kref) = (1 / (8.617e-5 * K)) - (1 / (8.617e-5 * Kref))
TtoK(T, Kref) =  1 / ((8.617e-5 * T) + (1 / Kref))

#boltzmann function
function boltz(B0::Float64,E::Float64,ΔT::Float64)
    B0 * exp(-E * ΔT)
end

#parameter structures
mutable struct Param 
    N::Int64
    r::Vector{Float64}
    a::Array{Float64,2}
    t_stop::Float64
end

#GLV derivatives
function dx!(dx,x,p,t)
    for i = 1:p.N
        dx[i] = x[i] * p.r[i]
        for j = 1:p.N
            dx[i] += x[i] * x[j] * p.a[i,j]           
        end 
    end
end


#Assembly callback functions

#1) condition - detect equilibrium
function detect_eq(integrator, abstol::Float64, reltol::Float64)
    x = !any(abs.(DiffEq.get_du(integrator)) .> abstol)
    y = (integrator.t - integrator.p.t_stop) > 1e3   
    return(x | y)
end

#2) affect! - Add new species
function add_species!(integrator, int_func::Function, grw_func::Function)
    #get arrays
    tmp = deepcopy(integrator.u)     

    #detect extinctions
    ext = tmp .> eps()
    #get new system size
    integrator.p.N = (1 + sum(ext))
    #resize integrator
    resize!(integrator,integrator.p.N)

    #put in old biomass
    integrator.u[1:(integrator.p.N-1)] .= tmp[ext]
    #put in new biomass
    integrator.u[end] = 0.01
    
    #update params
    #growth rate
    deleteat!(integrator.p.r, .!ext)
    push!(integrator.p.r, grw_func(1)[1])
    
    #interactions
    new_a = int_func(integrator.p.N, integrator.p.N)
    new_a[integrator.p.N, integrator.p.N] = -1.0
    #add old data  
    new_a[1:(integrator.p.N-1) , 1:(integrator.p.N-1) ] .= integrator.p.a[ext, ext]

    integrator.p.a = new_a
    integrator.p.t_stop = copy(integrator.t)
end

function add_at_equi(new_int::Function, grw_func::Function; abstol::Float64 = 1e-8, reltol::Float64 = 1e-6)
    condition = (u,t,integrator) -> detect_eq(integrator, abstol, reltol)
    affect! = (integrator) -> add_species!(integrator, new_int, grw_func)

    return(DiffEq.DiscreteCallback(condition, affect!; save_positions = (false, true)))
end

##predictions
function n_max(p)
    r_norm = p.r ./ mean(p.r)
    a_mean = sum(p.a,dims = 2)
end

function get_N(sol)
    length.(sol.u)
end

#gives probabiltiy of feasability as a function of richness
function p_feas_pred(N,ā,r̄,r_σ)
    #theta
    θ = (-ā*(N-1)) / (-ā*(N-1) + 1)
    
    #rnorm ditribution
    dr = LogNormal(-(r_σ^2)/2,r_σ)
    
    return((1 - cdf(dr,θ))^N)
end



#utility
function moving_avg(x::Vector, t_span::Int = 5)
    @assert t_span < length(x)
    
    avg = zeros(length(x))
    
    for i = eachindex(x)
        if i > t_span
            avg[i] = mean(x[(i-t_span):i])
        else
            avg[i] = mean(x[1:i])
        end
    end
    
    return(avg)
end

function mean_no_diag(x)
    a = sum(x,dims=2) .- [x[i,i] for i = size(x)[1]]
    b = mean(a) / (size(x)[1]-1)
    return(b)
end


function coef_var(x)
    std(x) / mean(x)
end
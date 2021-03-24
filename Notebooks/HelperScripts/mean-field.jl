
#mean-field
function dx_mean!(dx,x,p,t)
    for i = 1:p.N
        dx[i] = (x[i] * p.r[i]) + (x[i] * x[i] * p.a[i,i]) + (x[i] * mean(x) * p.N * p.a_bar)
    end
end


#for mean-field simulations
struct param_mean
    N::Int64
    r::Vector{Float64}
    a::Array{Float64,2}
    a_bar::Float64
end



function equi_app(p)
    a_ii = [a[i,i] for i = 1:p.N]
    r = p.r
    θ = mean(sum(p.a,dims = 2) - a_ii) ./ a_ii
    pred = (-r ./ a_ii) .+ (mean(r)./ a_ii) .* (θ ./ (θ .+ 1))
    return(pred)
end

function V(sol, t, meanfield)
    if meanfield
        equi_bio = equi_app(sol.prob.p)
    else
        equi_bio = sol[end]
    end
        V = mean(sol(t)) - mean(equi_bio) - mean(equi_bio) * log(mean(sol(t))/mean(equi_bio))
        return(V)
end   

function dV(sol,t,meanfield)
    if meanfield
        equi_bio = equi_app(sol.prob.p)
    else
        equi_bio = sol[end]
    end
    
    dV = sum((sol(t) .- equi_bio) .* (sol(t, Val{1}) ./ sol(t)))
    
    return(dV)
end
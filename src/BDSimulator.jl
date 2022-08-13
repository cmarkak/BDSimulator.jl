module BDSimulator

using Distributions
using DataFrames
using QuadGK
using Roots
using StaticArrays
using DifferentialEquations
using Random

export 
    simulateBDDeterministic,
    simulateBDStochastic

include("SSA.jl")

# Simulates a simple birth - death process of cells with b - a constant rate of birth and d*N a rate of 
# death linearly dependent on cell population using numberical solver for an ODE
function simulateBDDeterministic(;birth_rate::Float64=1.0,death_rate::Float64=0.1, timeSpan::Float64=100.0)
    f(n, p, t) = birth_rate - death_rate*n
    n_0 = 0
    tspan = (0.0,timeSpan)
    prob = ODEProblem(f,n_0,tspan)
    sol = solve(prob)
end

# Simulates a simple birth - death process of cells with b - a constant propensity of birth and d*N a propensity of 
# death linearly dependent on cell population using an implementation of the Gillespie SSA
function simulateBDStochastic(;birth_rate::Float64=1.0,death_rate::Float64=0.1, timeSpan::Float64=100.0)
    function F2(x,parms)
        (N,) = x
        (p_b,p_d) = parms
        [p_b,p_d*N]
    end

    x0 = [0]
    nu = reshape([[1];[-1]],2,1)
    parms = [birth_rate,death_rate]
    tf = timeSpan
    Random.seed!(1234)

    result = gillespie(x0,F2,nu,parms,tf)

    ssa_data(result)
end

end

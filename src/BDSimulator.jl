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
"
Simulates a simple birth - death process of cells with ba constant rate of birth - and d*N a rate of 
death linearly dependent on cell population using numberical solver for a the differential equation 
dN/dt = b - d*n

- **b** : the constant birth rate of the system (`Float64`).
- **d** : the death rate of the cells - dependent on the number of cells (`Float64`).
- **timeSpan** : the simulation time (`Float64`).
- **n0** : the initial number of cells (`Int64`).
"
function simulateBDDeterministic(;b::Float64=1.0,d::Float64=0.1, timeSpan::Float64=100.0, n0::Int64=1)
    f(n, p, t) = b - d*n
    n_0 = n0
    tspan = (0.0,timeSpan)
    prob = ODEProblem(f,n_0,tspan)
    sol = solve(prob)
end
"
Simulates a simple birth - death process of cells using the Gillespie stochastic simulation algorithm with propensities
p_b and and p_d for birth and death rates 

- **p_b** : the constant birth rate of the system (`Float64`).
- **p_d** : the death rate of the cells - dependent on the number of cells (`Float64`).
- **timeSpan** : the simulation time (`Float64`).
- **n0** : the initial number of cells (`Int64`).
"
function simulateBDStochastic(;p_b::Float64=1.0,p_d::Float64=0.1, timeSpan::Float64=100.0, n0::Int64=1)
    function F2(x,parms)
        (N,) = x
        (pb,pd) = parms
        [pb,pd*N]
    end

    x0 = [n0]
    nu = reshape([[1];[-1]],2,1)
    parms = [p_b,p_d]
    tf = timeSpan
    #Random.seed!(1234)

    result = gillespie(x0,F2,nu,parms,tf)

    ssa_data(result)
end

end

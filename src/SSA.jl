"A type storing the status at the end of a call to `ssa`:

- **termination_status** : whether the simulation stops at the final time (`finaltime`) or early due to zero propensity function (`zeroprop`)
- **nsteps** : the number of steps taken during the simulation.

"
struct SSAStats
    termination_status::String
    nsteps::Int64
end

"A type storing the call to `ssa`:

- **x0** : a `Vector` of `Int64`, representing the initial states of the system.
- **F** : a `Function` or a callable type, which itself takes two arguments; x, a `Vector` of `Int64` representing the states, and parms, a `Vector` of `Float64` representing the parameters of the system.
- **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : a `Vector` of `Float64` representing the parameters of the system.
- **tf** : the final simulation time (`Float64`).
- **alg** : the algorithm used (`Symbol`, either `:gillespie`, `jensen`, or `tjc`).
- **tvc** : whether rates are time varying.
"
struct SSAArgs{X,Ftype,N,P}
    x0::X
    F::Ftype
    nu::N
    parms::P
    tf::Float64
    alg::Symbol
    tvc::Bool
end

"
This type stores the output of `ssa`, and comprises of:

- **time** : a `Vector` of `Float64`, containing the times of simulated events.
- **data** : a `Matrix` of `Int64`, containing the simulated states.
- **stats** : an instance of `SSAStats`.
- **args** : arguments passed to `ssa`.

"
struct SSAResult
    time::Vector{Float64}
    data::Matrix{Int64}
    stats::SSAStats
    args::SSAArgs
end


"
This function is a substitute for `StatsBase.sample(wv::WeightVec)`, which avoids recomputing the sum and size of the weight vector, as well as a type conversion of the propensity vector. It takes the following arguments:

- **w** : an `Array{Float64,1}`, representing propensity function weights.
- **s** : the sum of `w`.
- **n** : the length of `w`.

"
function pfsample(w::AbstractArray{Float64,1},s::Float64,n::Int64)
    t = rand() * s
    i = 1
    cw = w[1]
    while cw < t && i < n
        i += 1
        @inbounds cw += w[i]
    end
    return i
end

"
This function performs Gillespie's stochastic simulation algorithm. It takes the following arguments:

- **x0** : a `Vector` of `Int64`, representing the initial states of the system.
- **F** : a `Function` or a callable type, which itself takes two arguments; x, a `Vector` of `Int64` representing the states, and parms, a `Vector` of `Float64` representing the parameters of the system.
- **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : a `Vector` of `Float64` representing the parameters of the system.
- **tf** : the final simulation time (`Float64`).
"
function gillespie(x0::AbstractVector{Int64},F::Base.Callable,nu::AbstractMatrix{Int64},parms::AbstractVector{Float64},tf::Float64)
    # Args
    args = SSAArgs(x0,F,nu,parms,tf,:gillespie,false)
    # Set up time array
    ta = Vector{Float64}()
    t = 0.0
    push!(ta,t)
    # Set up initial x
    nstates = length(x0)
    x = copy(x0')
    xa = copy(Array(x0))
    # Number of propensity functions
    numpf = size(nu,1)
    # Main loop
    termination_status = "finaltime"
    nsteps = 0
    while t <= tf
        pf = F(x,parms)
        # Update time
        sumpf = sum(pf)
        if sumpf == 0.0
            termination_status = "zeroprop"
            break
        end
        dt = rand(Exponential(1/sumpf))
        t += dt
        push!(ta,t)
        # Update event
        ev = pfsample(pf,sumpf,numpf)
        if x isa SVector
            @inbounds x[1] += nu[ev,:]
        else
            deltax = view(nu,ev,:)
            for i in 1:nstates
                @inbounds x[1,i] += deltax[i]
            end
        end
        for xx in x
            push!(xa,xx)
        end
        # update nsteps
        nsteps += 1
    end
    stats = SSAStats(termination_status,nsteps)
    xar = transpose(reshape(xa,length(x),nsteps+1))
    return SSAResult(ta,xar,stats,args)
end


"This takes a single argument of type `SSAResult` and returns a `DataFrame`."
function ssa_data(s::SSAResult)
    hcat(DataFrame(time=s.time),DataFrame(s.data,:auto))
end

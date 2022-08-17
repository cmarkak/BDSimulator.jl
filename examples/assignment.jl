push!(LOAD_PATH,"../src/")
#using BDSimulator

sweep_range = range(-2,2)

deterministic_sims = dict()
stochastic_sims = dict()

for expon in sweep_range
    stochastic_sims["b=$(expon)"] = simulateBDStochastic(b=10.0^expon, timeSpan=1000.0)
    stochastic_sims["a=$(expon)"] = simulateBDStochastic(a=10.0^expon, timeSpan=1000.0)
    deterministic_sims["b=$(expon)"] = simulateBDDeterministic(b=10.0^expon, timeSpan=1000.0)
    deterministic_sims["a=$(expon)"] = simulateBDDeterministic(a=10.0^expon, timeSpan=1000.0)
end
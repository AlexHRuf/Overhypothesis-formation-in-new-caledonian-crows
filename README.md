# Overhypothesis formation in new caledonian crows
Modeling the formation of overhypotheses in New Caledonian Crows. Crows sample different food types (high- and low-value) and learn about their distribution (represented by α and β) within and across a series of containers. In the test condition (final container), this abstract knowledge informs the probability of finding a high-value item (θ[high-value]).

![me](https://github.com/AlexHRuf/Overhypothesis-formation-in-new-caledonian-crows/blob/main/Animations/theta_evolution_combined_dark.gif)
Mixed group (left): Crows in the mixed group learn that there is an equal amount of low- and high-value food across containers (β = 0.5) and that this distribution is reflected in each container (high α).
Uniform group (right): Crows in the uniform group learn that there is an equal amount of low- and high-value food across containers (β = 0.5) but that only one item type can be found in each container (low α).


# How to run the .ipynb scripts using Julia
First, you will need to install Julia: **https://julialang.org/downloads/**
Once julia is installed on your machine, you can run it by entering "julia" in the terminal. For further information about getting started in Julia and how it differes from other languages see: **https://docs.julialang.org/en/v1/manual/getting-started/**. 

Next, you will need to install IJulia by following the instructions here: **https://github.com/JuliaLang/IJulia.jl**.
Once successfully installed you will be able to run a Julia kernel in Jupyter Notebook.

To run the scripts in this repository first install the dependencies by running:
```julia
using Pkg

dependencies = ["Turing", "Distributions", "Random", "MCMCChains", "Plots", "StatsPlots", "Measures", "BSON"]

Pkg.add(dependencies)
```
We are using Turing.jl to define hierarchical bayesian models. For an introduction to Turing.jl see: **https://turinglang.org/docs/tutorials/00-introduction/**.

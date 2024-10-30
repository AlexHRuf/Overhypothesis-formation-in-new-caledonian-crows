# Overhypothesis formation in new caledonian crows
Modelling the formation of overhypotheses in New Caledonian Crows. Crows [in this example the ideal model of a crow capable of abstraction] sample different food types (high- and low-value) and learn about their distribution within and across a series of containers (represented by α and β). In the test condition (final container, shown below), this abstract knowledge informs the probability of finding a high-value item (θ[high-value]), while consecutively sampling low-value items.

![me](https://github.com/AlexHRuf/Overhypothesis-formation-in-new-caledonian-crows/blob/main/Models/Animations/theta_evolution_groups_tests.gif)

When continuously sampling low-value items from the Test 1 container (red), mixed group birds predict a high probability of finding a high-value item even after several samples. Meanwhile, uniform group birds only need to see one low-value sample to be almost certain that they will not find a high-value item in Test 1. The probability of sampling a high-value item from the Test 2 container (gray) stays relatively constant. This creates the incentive for crows to switch from the Test 1 to the Test 2 container. Crows that are in the uniform group are expected to make the switch earlier (after fewer samples) than crows from the mixed group.

# How to run the .ipynb scripts using Julia
First, you will need to install Julia: **https://julialang.org/downloads/**.
Once Julia is installed on your machine, you can run it by entering "julia" in the terminal. For further information on getting started in Julia and how it differs from other languages see: **https://docs.julialang.org/en/v1/manual/getting-started/**. 

Next, you will need to install IJulia by following the instructions here: **https://github.com/JuliaLang/IJulia.jl**.
Once successfully installed you will be able to run a Julia kernel in Jupyter Notebook.

To run the scripts in this repository first install the dependencies by running:
```julia
using Pkg

dependencies = ["Turing", "Distributions", "Random", "MCMCChains", "Plots", "StatsPlots", "Measures", "BSON"]

Pkg.add(dependencies)
```
We are using Turing.jl to define hierarchical bayesian models. For an introduction to Turing.jl see: **https://turinglang.org/docs/tutorials/00-introduction/**.

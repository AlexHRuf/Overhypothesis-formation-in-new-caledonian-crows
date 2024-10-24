# Overhypothesis formation in new caledonian crows
Modeling the formation of overhypotheses in a hidden food sampling task.

![me](https://github.com/AlexHRuf/Overhypothesis-formation-in-new-caledonian-crows/blob/main/Animations/theta_evolution_test1_mixed.gif)

# How to run the model scripts
First, you will need to install Julia: **https://julialang.org/downloads/**
Once julia is installed on your machine, you can run it by enterin "julia" in the terminal.

Next, you will need to install IJulia by following the instructions here: **https://github.com/JuliaLang/IJulia.jl**
Once successfully installed you will be able to run Julia in Jupyter Notebook.

To run the scripts in this repository first install the dependencies by running:
```julia
using Pkg

dependencies = ["Turing", "Distributions", "Random", "MCMCChains", "Plots", "StatsPlots", "Measures", "BSON"]

Pkg.add(dependencies)
```

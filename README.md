# Forced Randomization

Julia code for simulating four different randomization schemes in multicenter 2-arm clinical trials with drug supply in mind.

## Overview:

- funcs.jl contains essential functions necessary to run main.jl, like creating the patient and treatment lists, and assigning scenario values. It also contains a few parameters which can be altered relating to the resupply strategies.
- Simulation.jl contains the logic behind each scenario and is where each specification is actually carried out.
- plotting.jl contains a few functions used to plot and export the results.
- main.jl is the main body for the code. Notably, it contains the constant simulation parameters which can be altered to test different sample sizes, numbers of centers, and recruitment rates. Additionally, it is where all the data is gathered, stored, and exported.

- Running include("main.jl") in a Julia REPL will output results described in section "An Example"

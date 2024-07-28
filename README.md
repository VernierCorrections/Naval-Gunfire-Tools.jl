This is a revision of my earlier external ballistics script, rewritten entirely in the Julia programming language, and utilizing a Vern9 ODE solver. Admittedly this version has room for improvement via more extensive use of vectorized code, but should run faster than my previous MATLAB script, and possesses new features (such as displaying a range table).

# To-Do
Higher-order gravity terms using automatic differentiation.

# Sources
In addition to any sources cited on my previous project, I utilized Vincenty's 1975 paper, in addition to SciML's DifferentialEquations.jl package, in addition to the Plots.jl, AstroTime.jl, and Interpolations.jl packages.

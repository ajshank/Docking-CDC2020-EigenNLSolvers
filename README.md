# Docking-CDC2020-EigenNLSolvers
Independent code to solve the non-linear equations for trajectory generation
derived in:
> "Dynamic Path Generation for Multirotor Aerial Docking in Forward Flight",   
> A. Shankar, S. Elbaum, and C. Detweiler;   
> Conference on Decision & Control (CDC), 2020.

The problem has 6 non-linear algebraic equations with 6 unknowns. Two classic
solvers/algorithms are used here from the Eigen `unsupported` library: 
Trust-Region, and Levenberg-Marquardt. Both produce identical results.

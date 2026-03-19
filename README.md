# Shallow Water Equations

## Project Scope


- Derive the two-dimensional linearized shallow water equations (SWE) as mathematical model.
    
- Derive and implement a suitable numerical solver in C++.  
A discontinuous Galerkin (DG) method is a strong option, since it is well suited for hyperbolic PDEs and parallelization and I have experience implementing it for Poisson Equ. The project is not restricted to DG if another method turns out to be more appropriate.
    
- Verify correctness of the implementation against known *analytical* solutions for the linearized problem, including *convergence* studies.
    
- Parallelize the solver in order to study performance and *scalability*, for example using `OpenMP` or `CUDA`.

- Extend the model to the full two-dimensional __nonlinear__ shallow water equations as a natural continuation once the linearized solver is established.

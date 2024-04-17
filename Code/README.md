This folder contains the two folders for the code and results of the porous medium equation 

```math
\partial_t u = \Delta(u^m) + \beta u,
```

and the biofilm model

```math
\begin{cases}
    \partial_t u &= \Delta(\Phi(u)) + f(v)u, \\
    \partial_t v &= \nabla \cdot \left(D(u)\nabla v\right) + g(u,v).
\end{cases}
```
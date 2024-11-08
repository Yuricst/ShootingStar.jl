# ShootingStar.jl

This package implements the two-stage shooting method for trajectory design in Newtonian systems with impulsive control. 

The trajectory is computed by discretizing it into `N` nodes, composed of `N-1` segments. 
The algorithm consists of an inner-loop that drives the position discrepancies at each node to 0 by modifying the velocities, and an outer-loop that minimizes the velocity discrepancies by modifying the positions of the intermediate nodes. 

Note that the inner-loop has `3(N-1)` degrees of freedom for `3(N-1)` constraints, resulting in square Newton-Raphson process, while the outer-loop has `3(N-2)` degrees of freedom and `3N` constraints, resulting in a least-squares update. 
Due to the least-squares nature of the outer-loop, the resulting trajectory is minimum-energy (as opposed to minimum-fuel). 

## Quick start

1. `git clone` this repository and `cd` to its root
2. Start julia REPL and activate environment 

```julia-repl
pkg> activate .
```

3. Run tests to make sure installation is successfull

```julia-repl
(ShootingStar) pkg> test
```


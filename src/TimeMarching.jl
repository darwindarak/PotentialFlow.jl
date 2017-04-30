"""
A collection of time marching algorithms for vortex modeling

## Exports

- Forward Euler (`forward_euler!`)
- 4th-order Runge-Kutta (`rk4!`)

All schemes have the following function signature

    scheme_name!(x₊, x, t, Δt, compute_ẋ!, update_x!, k)

where

- `x₊`: where the next state is stored
- `x`: current state
- `t`: current time
- `Δt`: time step
- `compute_ẋ!(k, x, t)`: a function with computes the velocity of the state `x` at time `t` and stores the result in `k`
- `update_x!(x₊, x, Δt, k...)`: a function that updates the current state `x` with velocity `k` over the time interval `Δt` and stores the result in `x₊`
- `k...`: one or more temporary variables to store the velocity (see the docstrings for specific schemes)
"""
module TimeMarching

export forward_euler!, rk4!

"""
    forward_euler!(x₊, x, t, Δt, compute_ẋ!, update_x!, k)

Update the state `x` at time `t` using forward Euler

Expects only a single `k`, representing the velocity at time `t`.
"""
function forward_euler!(x₊, x, t, Δt, compute_ẋ!, update_x!, k)
    compute_ẋ!(k, x, t)
    update_x!(x₊, x, Δt, k)
    return x₊, Δt
end

"""
    rk4!(x₊, x, t, Δt, compute_ẋ!, update_x!, k)

Update the state `x` at time `t` using 4th-order Runge-Kutta

Expects `k` to be indexable from 1 to 4, representing velocities at the Runge-Kutta half-steps.
"""
function rk4!(x₊, x, t, Δt, compute_ẋ!, update_x!, k)
    compute_ẋ!(k[1], x, t)

    update_x!(x₊, x, 0.5Δt, k[1])
    compute_ẋ!(k[2], x₊, t + 0.5Δt)

    update_x!(x₊, x, 0.5Δt, k[2])
    compute_ẋ!(k[3], x₊, t + 0.5Δt)

    update_x!(x₊, x, Δt, k[3])
    compute_ẋ!(k[4], x₊, t + Δt)

    ẋ = @. (k[1] + 2k[2] + 2k[3] + k[4])/6

    update_x!(x₊, x,  Δt, ẋ)
    return x₊, Δt
end

end

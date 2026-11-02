# Linear algebra

Expressions implement [VectorInterface.jl](https://github.com/Jutho/VectorInterface.jl), with the inner product `inner(x, y) = scalar(x' * y)`. This requires `adjoint` to be defined for your types, and `x' * y` must sort to an expression with no factors left (apart from pairs of terms skipped by keys, see [Fast inner products](@ref)), so that [`scalar`](@ref NonCommutativeProducts.scalar) can extract its value. `inner` sorts the pairs it doesn't skip, whether or not autosort is enabled.

With [KrylovKit.jl](https://github.com/Jutho/KrylovKit.jl) loaded, expressions can also be used as linear operators acting by left multiplication. These features live in a package extension, which `using KrylovKit` loads.

Example: Let's define kets and bras and exponentiate a Hamiltonian:
```jldoctest kets
julia> using KrylovKit, LinearAlgebra

julia> import NonCommutativeProducts as NC

julia> NC.enable_autosort!();

julia> struct Ket
           n::Int
           bra::Bool
           Ket(n::Int, bra::Bool=false) = new(n, bra)
       end

julia> Base.show(io::IO, k::Ket) = print(io, k.bra ? "⟨" : "|", k.n, k.bra ? "|" : "⟩")

julia> Base.adjoint(k::Ket) = Ket(k.n, !k.bra)

julia> NC.@nc Ket

julia> NC.mul_effect(a::Ket, b::Ket) = a.bra && !b.bra ? Int(a.n == b.n) : nothing

julia> k0, k1 = Ket(0), Ket(1);

julia> H = k0 * k1' + k1 * k0'
|0⟩*⟨1| + |1⟩*⟨0|

julia> vals, vecs = eigsolve(H, k0, 2, :SR; ishermitian=true);

julia> x, info = exponentiate(1im*H, pi/2, k0);

julia> norm(x - im*k1) < 1e-12
true
```

## Fast inner products

`inner(x, y)` forms `a' * b` for every pair of a term `a` of `x` and a term `b` of `y`, so its cost is proportional to the product of the numbers of terms. With [`contraction_key`](@ref NonCommutativeProducts.contraction_key) defined, `inner` skips the pairs whose keys show that they are orthogonal, which takes time proportional to the number of terms when each term matches few others.

Keys define orthogonality: `inner` skips pairs whose keys disagree, whatever `mul_effect` gives and without sorting them first, so terms must list subsystems in the same order. Here `⟨n|m⟩` is zero unless `n == m`, so kets are keyed by their label:
```jldoctest kets
julia> NC.contraction_key(k::Ket) = k.bra ? nothing : k.n

julia> using NonCommutativeProducts.VectorInterface: inner

julia> inner(k0 + 2k1, k1)
2
```

## Truncation with an error bound

[`truncate_terms`](@ref NonCommutativeProducts.truncate_terms) drops small terms of an expression and returns a bound on the error, valid in every representation where the atom bounds hold. For `x = Σᵢ cᵢ Pᵢ`, dropping a set `D` of terms gives `x̃` with
```
‖x - x̃‖ ≤ Σ_{i ∈ D} |cᵢ| b(Pᵢ)
```
by the triangle inequality, where `b(P)` is the product of the operator norm bounds of the factors of `P`. The true error can be much smaller.

To make terms droppable, define [`opnorm_bound`](@ref NonCommutativeProducts.opnorm_bound) for your atoms. Terms containing an atom without a bound, such as a boson, are kept.

`rtol` is relative to the total weight of the droppable terms, kept or not. That total can far exceed `‖x‖` and leaves out terms without a bound, so `rtol` is not a relative error in `‖x‖`. Combine like terms before truncating, e.g. with autosort, since terms that would cancel otherwise use up the budget.

For the kets above, normalized kets and bras have norm 1:
```jldoctest kets
julia> NC.opnorm_bound(::Ket) = 1

julia> y, err = NC.truncate_terms(H + 0.01 * k0 * k0' + 0.001 * k1 * k1'; atol=0.005);

julia> y
0.01*|0⟩*⟨0| + |0⟩*⟨1| + |1⟩*⟨0|

julia> err
0.001
```

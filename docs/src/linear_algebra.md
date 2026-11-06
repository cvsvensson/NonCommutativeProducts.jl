# Linear algebra

Expressions implement [VectorInterface.jl](https://github.com/Jutho/VectorInterface.jl), with the inner product `inner(x, y) = scalar(x' * y)`. This requires `adjoint` to be defined for your types, and `x' * y` must sort to an expression with no factors left (apart from pairs of terms whose keys differ, see [Fast inner products](@ref)), so that [`scalar`](@ref NonCommutativeProducts.scalar) can extract its value. `inner` sorts the pairs it doesn't skip, also inside `with_autosort(false)`. It sorts `x` and `y` too unless autosort is on and they can have keys, so `mul_effect` must be defined for their neighbouring factors.

With [KrylovKit.jl](https://github.com/Jutho/KrylovKit.jl) loaded, expressions can also be used as linear operators acting by left multiplication. These features live in a package extension, which `using KrylovKit` loads. Operators given as plain functions should return sorted results, because unmerged terms can make solvers report wrong eigenvalues as converged; `x -> sort(H * x)` does this when autosort is off.

Example: Let's define kets and bras and exponentiate a Hamiltonian:
```jldoctest kets
julia> using KrylovKit, LinearAlgebra

julia> import NonCommutativeProducts as NC

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

Terms that first differ at kets of one subsystem are orthogonal; otherwise `inner` computes the full product, so any order is correct and one order is fastest. Mixing orders can be tens of times slower.

`mul_effect(k', l)` decides whether kets are in one subsystem. See [`contraction_key`](@ref NonCommutativeProducts.contraction_key) for what it must return.

Sort a sum formed inside `with_autosort(false)` before calling `inner` with autosort on: equal products left in separate terms cancel only up to roundoff, so a zero vector gets a norm around `1e-8`.

Here `⟨n|m⟩` is zero unless `n == m`, so kets are keyed by their label:
```jldoctest kets
julia> NC.contraction_key(k::Ket) = k.bra ? nothing : k.n

julia> using NonCommutativeProducts.VectorInterface: inner

julia> inner(k0 + 2k1, k1)
2
```

## Truncation with an error bound

[`truncate_terms`](@ref NonCommutativeProducts.truncate_terms) drops small terms of an expression and returns a bound on the error. To make terms droppable, define [`opnorm_bound`](@ref NonCommutativeProducts.opnorm_bound) for your atoms. Sort an expression formed inside `with_autosort(false)` before truncating, since terms that would cancel otherwise use up the budget.

For the kets above, normalized kets and bras have norm 1:
```jldoctest kets
julia> NC.opnorm_bound(::Ket) = 1

julia> y, err = NC.truncate_terms(H + 0.01 * k0 * k0' + 0.001 * k1 * k1'; atol=0.005);

julia> y
0.01*|0⟩*⟨0| + |0⟩*⟨1| + |1⟩*⟨0|

julia> err
0.001
```

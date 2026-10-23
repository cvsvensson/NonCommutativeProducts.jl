# Linear algebra

Expressions implement [VectorInterface.jl](https://github.com/Jutho/VectorInterface.jl), with the inner product `inner(x, y) = scalar(x' * y)`. This requires `adjoint` to be defined for your types, and `x' * y` must sort to an expression with no factors left, so that [`scalar`](@ref NonCommutativeProducts.scalar) can extract its value. `inner` sorts `x' * y` whether or not autosort is enabled.

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

`inner(x, y)` forms `a' * b` for every pair of a term `a` of `x` and a term `b` of `y`, so its cost is proportional to the product of the numbers of terms. Defining [`contraction_key`](@ref NonCommutativeProducts.contraction_key) for your atom types makes `inner` form only the products of terms whose keys are equal, which takes time proportional to the number of terms when each key matches few terms. The keys must be equal whenever `a' * b` can be nonzero. A term with a factor whose key is `nothing`, the default, is paired with every term.

For the kets above, `⟨n|m⟩` is zero unless `n == m`, so the label is a key:
```jldoctest kets
julia> NC.contraction_key(k::Ket) = k.n

julia> using NonCommutativeProducts.VectorInterface: inner

julia> inner(k0 + 2k1, k1)
2
```

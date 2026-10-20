# Traversing expressions

```@meta
DocTestSetup = quote
    using NonCommutativeProducts
    import NonCommutativeProducts: @nc, mul_effect
    struct Fermion
        label::Int
        dagger::Bool
    end
    Base.adjoint(x::Fermion) = Fermion(x.label, !x.dagger)
    Fermion(k) = Fermion(k, false)
    Base.show(io::IO, x::Fermion) = print(io, "c", x.dagger ? "†" : "", "[", x.label, "]")
    @nc Fermion
    function mul_effect(a::Fermion, b::Fermion)
        (a.dagger, a.label) == (b.dagger, b.label) && return 0
        (!a.dagger, a.label) < (!b.dagger, b.label) && return nothing
        (!a.dagger, a.label) == (b.dagger, b.label) && return -b*a + 1
        return -b*a
    end
    NonCommutativeProducts.disable_autosort!()
end
```

The examples on this page use the `Fermion` type defined on the [Home](@ref "Example: Fermions") page.

[`ncmapreduce(f, (add, mul), x; scalarmap)`](@ref NonCommutativeProducts.ncmapreduce) maps `f` over every factor in `x`, combines the factors in each term with `mul`, and then combines the terms with `add`. Scalars are passed through `scalarmap`. For example, the highest number of factors in any term of a fermion expression is
```jldoctest traversal
julia> import NonCommutativeProducts: ncmapreduce

julia> op = 1 + 2Fermion(1)'*Fermion(2) + Fermion(1)'*Fermion(2)'*Fermion(3)
1I + 2*c†[1]*c[2] + c†[1]*c†[2]*c[3]

julia> ncmapreduce(x -> 1, (max, +), op; scalarmap=zero)
3
```
With the ordinary `(+, *)`, one can substitute a matrix for each factor:
```jldoctest traversal
julia> using LinearAlgebra

julia> c = [0 1; 0 0];

julia> ncmapreduce(x -> x.dagger ? c' : c, (+, *), Fermion(1)'*Fermion(1) + 2; scalarmap=x -> x*I)
2×2 Matrix{Int64}:
 2  0
 0  3
```

```@meta
DocTestSetup = nothing
```

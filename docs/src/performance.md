# Performance tips

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
    NonCommutativeProducts.disable_autosort!()
end
```

The examples on this page use the `Fermion` type defined on the [Home](@ref "Example: Fermions") page.

This package is flexible, but not very efficient. Sorting is done via bubble sort, which is convenient for this use case since it is based on repeatedly swapping adjacent elements where commutation relations can be used. But it does not scale well with the length of the list, so it won't perform well for products of many elements.

```julia
# Example timing: 108.298 ms (1142139 allocations: 78.37 MiB)
op = prod(Fermion(n) + Fermion(n)' + 1 for n in 1:10)
#= Sum with 59049 terms:
1I-c†[10]*c[1]-c†[10]*c[1]*c[2]*c[3]-c†[10]*c[1]*c[2]*c[3]*c[4]*c[5] + ...
=#
```

To cut down on allocations, one can use [`NonCommutativeProducts.add!!`](@ref) which tries to perform addition in place, but widens the type if not possible.
```julia
op = 0
for n in 1:100
    op += Fermion(n)'*Fermion(n)
end
# Example timing: 105.000 μs (2221 allocations: 543.15 KiB)

op2 = zero(op)
for n in 1:100
    op2 = NonCommutativeProducts.add!!(op2, Fermion(n)'*Fermion(n))
end
# Example timing: 44.300 μs (1613 allocations: 108.64 KiB)
op == op2 # true
```

Mixing different symbolic types together is slower than using a single type.

Disable autosort to improve performance if you don't want automatic sorting.

Finally, instead of using ordinary addition and multiplication in `mul_effect`, you can use [`Swap`](@ref NonCommutativeProducts.Swap) and [`AddTerms`](@ref NonCommutativeProducts.AddTerms) which are more efficient. Their behaviour is:
* `Swap(λ::Number)`: Replaces `a*b` by `λ*b*a`.
* `AddTerms(terms)`: `a*b` should be replaced by a sum of terms. `terms` should be an iterable such as a vector or a tuple, and the elements can be of the other allowed return types except `nothing`.

For example, the `mul_effect` for fermions can be implemented as
```jldoctest performance; output = false, filter = r"\d+ methods?"
import NonCommutativeProducts: Swap, AddTerms
function mul_effect(a::Fermion, b::Fermion)
    (a.dagger, a.label) == (b.dagger, b.label) && return 0
    (!a.dagger, a.label) < (!b.dagger, b.label) && return nothing
    (!a.dagger, a.label) == (b.dagger, b.label) && return AddTerms((Swap(-1), 1))
    return Swap(-1)
end

# output

mul_effect (generic function with 1 method)
```
which sorts the same way:
```jldoctest performance
julia> Fermion(1)*Fermion(1)' |> sort
1I - c†[1]*c[1]
```

```@meta
DocTestSetup = nothing
```

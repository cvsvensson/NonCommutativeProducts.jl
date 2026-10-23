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

Each product is sorted with bubble sort, which repeatedly swaps adjacent factors so that commutation relations apply. The time to sort one product grows quadratically with its number of factors, so long products are slow.

A product of sums is formed one pair of terms at a time: each pair is sorted on its own and added straight into the result. The time is proportional to the number of term pairs and the memory to the size of the result. A product that sorts to zero allocates nothing, and each term added to the result costs one copy of its factors.

```julia
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

op2 = zero(op)
for n in 1:100
    op2 = NonCommutativeProducts.add!!(op2, Fermion(n)'*Fermion(n))
end
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

When products of a term with a sum have structure that `mul_effect` cannot use, such as an operator string that acts on every basis state of a sum of kets with one precomputation, define a method of [`mul_terms!!`](@ref NonCommutativeProducts.mul_terms!!) for your types. With autosort enabled, `*` calls `mul_terms!!(c, x, b)` for each term `x` of the left factor, with autosort disabled during the call. A method reads the terms of `b` with [`termpairs`](@ref NonCommutativeProducts.termpairs) and adds the results to `c` with [`add!!`](@ref NonCommutativeProducts.add!!).

For inner products of large sums, define [`contraction_key`](@ref NonCommutativeProducts.contraction_key) for atoms that behave as basis vectors, such as kets keyed by their basis label, so that `inner` forms `a' * b` only for pairs of terms whose factor keys are equal; see [Fast inner products](@ref).

```@meta
DocTestSetup = nothing
```

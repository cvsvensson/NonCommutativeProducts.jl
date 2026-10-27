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

When products of a term with a sum have structure that `mul_effect` cannot use, define a method of [`mul_terms!!`](@ref NonCommutativeProducts.mul_terms!!) for your types. With autosort enabled, `mul!!` calls `mul_terms!!(c, x, b)` for each term `x` of the left factor, and so does `*` when either factor is a sum. Autosort is disabled during the call, so the method must add terms that are already sorted.

For example, `|n⟩⟨m|` applied to a sum of kets is the coefficient of `|m⟩` times `|n⟩`. A method can find it in one pass over the terms of the sum with [`termpairs`](@ref NonCommutativeProducts.termpairs), and add it to `c` with [`add!!`](@ref NonCommutativeProducts.add!!). Other products go through the default method:
```jldoctest multerms
julia> import NonCommutativeProducts as NC

julia> NC.enable_autosort!();

julia> struct Ket
           n::Int
           bra::Bool
       end

julia> Base.adjoint(k::Ket) = Ket(k.n, !k.bra)

julia> Base.show(io::IO, k::Ket) = print(io, k.bra ? "⟨" : "|", k.n, k.bra ? "|" : "⟩")

julia> NC.@nc Ket

julia> NC.mul_effect(a::Ket, b::Ket) = a.bra && !b.bra ? Int(a.n == b.n) : nothing

julia> isket(fs) = length(fs) == 1 && !only(fs).bra;

julia> function NC.mul_terms!!(c, x::NC.NCMul{<:Any,Ket}, b::NC.NCAdd)
           fs = NC.factors(x)
           if length(fs) != 2 || fs[1].bra || !fs[2].bra || !all(isket ∘ first, NC.termpairs(b))
               return invoke(NC.mul_terms!!, Tuple{Any,NC.NCMul,NC.NCAdd}, c, x, b)
           end
           for (kets, coeff) in NC.termpairs(b)
               if only(kets).n == fs[2].n
                   c = NC.add!!(c, NC.coefficient(x) * coeff * Ket(fs[1].n, false))
               end
           end
           return c
       end;

julia> ψ = sum(n * Ket(n, false) for n in 1:3);

julia> Ket(1, false) * Ket(3, true) * ψ
3*|1⟩
```

For inner products of large sums, see [Fast inner products](@ref).

```@meta
DocTestSetup = nothing
```

# Grades are metadata and never affect the value. A nonzero expression has the grades of its terms, while a zero has
# the grades of how it was computed (see `grades`). A zero sum stores these in `zero_grades`, with `nothing` for ∅.
# Invariant: every NCAdd returned or mutated by a public operation has correct `zero_grades` if it is zero, and
# `zero_grades === nothing` otherwise. Public operations restore it at exit with _with_zero_grades; internal helpers
# don't maintain it. Grade sets are Set{Any} and never mutated once stored.
#
# The grades of a nonzero destination of add!!/mul!! are ignored when the result is zero. This is sound if mul_effect
# respects the grading (see check_grading): c + v == 0 with c != 0 implies that the grades of c are among those of v.

"""
    grade(atom)

The grade of `atom`. Define it for your types to make them graded; the default is `nothing`, the identity grade.
Every term of [`mul_effect`](@ref)`(a, b)` must have the product of the grades of `a` and `b`, see
[`check_grading`](@ref). Within one grading, this is [`mul_grade`](@ref)`(grade(a), grade(b))`. Grades of different
types are combined according to their [`grading`](@ref).

For other expressions, `grade(x) == only(grades(x))`, and it throws an `ArgumentError` if there are several grades or
none.

# Example
```julia
# Fermions with charge as the grade
NonCommutativeProducts.grade(f::Fermion) = f.creation ? 1 : -1
NonCommutativeProducts.mul_grade(g::Int, h::Int) = iszero(g + h) ? nothing : g + h
NonCommutativeProducts.adjoint_grade(g::Int) = -g
```
"""
grade(x) = nothing
grade(x::Number) = iszero(x) ? throw(_no_grade_error(x, Set{Any}())) : nothing
grade(x::UniformScaling) = grade(x.λ)
function grade(x::NCMul{C,S}) where {C,S}
    isempty(x.factors) && return grade(prefactor(x))
    _is_ungraded(S) && return nothing
    return foldl((g, f) -> _mul_grade(g, grade(f)), x.factors; init=nothing)
end
function grade(x::NCAdd)
    gs = _grades(x)
    length(gs) == 1 || throw(_no_grade_error(x, gs))
    return only(gs)
end

function _no_grade_error(x, gs)
    isempty(gs) && return ArgumentError("$x is a zero without grades (e.g. the number 0 or a zero of unknown origin), so it has no grade")
    return ArgumentError("$x has several grades $gs; use `grades`")
end

"""
    mul_grade(g, h)

The product of the grades `g` and `h` of one grading, neither of which is `nothing`. It must be associative, and
return `nothing` for the identity grade or else a grade of the same grading.

It is only called for grades of one grading: grades of the same type, or of different types with the same
[`grading`](@ref). So a family of types sharing a key needs `mul_grade` methods between its member types. Grades of
different gradings are combined into a [`ProductGrade`](@ref) instead.
"""
function mul_grade end

"""
    adjoint_grade(g)

The grade of the adjoint of an expression of grade `g` (not `nothing`). It must return `nothing` or a grade of the same
grading. Only needed for adjoints of zero sums.
"""
function adjoint_grade end

"""
    grading(G::Type)

The key of the grading that grades of type `G` belong to, or `nothing` if it is undeclared (the default). Keys are
compared with `isequal`, and should be types or symbols.

Grades of the same type are multiplied with [`mul_grade`](@ref), and need no declaration. Grades of different types
are combined by their keys:
- equal keys: one grading, multiplied with `mul_grade`, which must be defined between the types;
- different keys: independent gradings, combined into a [`ProductGrade`](@ref);
- an undeclared key: an `ArgumentError`.

# Example
```julia
NonCommutativeProducts.grading(::Type{Charge}) = Charge
NonCommutativeProducts.grading(::Type{<:U1}) = U1 # one grading spread over the types U1{T}
```
"""
grading(::Type) = nothing

_mul_grade(::Nothing, ::Nothing) = nothing
_mul_grade(::Nothing, h) = h
_mul_grade(g, ::Nothing) = g
function _mul_grade(g, h)
    (g isa ProductGrade || h isa ProductGrade) && return _mul_product_grade(g, h)
    typeof(g) === typeof(h) && return mul_grade(g, h)
    kg, kh = grading(typeof(g)), grading(typeof(h))
    (isnothing(kg) || isnothing(kh)) && throw(_undeclared_grading_error(g, h))
    isequal(kg, kh) && return mul_grade(g, h)
    return ProductGrade((g, h), Val(:unchecked))
end
_adjoint_grade(::Nothing) = nothing
_adjoint_grade(g) = adjoint_grade(g)

@noinline function _undeclared_grading_error(g, h)
    Ts = unique!([typeof(x) for x in (g, h) if !(x isa ProductGrade) && isnothing(grading(typeof(x)))])
    defs = join(("`NonCommutativeProducts.grading(::Type{$T})`" for T in Ts), " and ")
    return ArgumentError("cannot multiply the grades $(repr(g)) and $(repr(h)) of different types: define $defs to \
        treat them as independent gradings, or, if they are the same grading, give them the same key and define \
        `mul_grade` between them")
end

# Whether `c`, the result of mul_grade or adjoint_grade within the grading `key`, stays in that grading
_is_closed(c, key) = isnothing(c) || isequal(grading(typeof(c)), key)

@noinline function _closure_error(f, args, c)
    call = "$(nameof(f))($(join(map(repr, args), ", ")))"
    return ArgumentError("$call returned $(repr(c)), which is neither `nothing` nor a grade of the grading \
        $(repr(grading(typeof(first(args)))))")
end

"""
    ProductGrade(g1, g2, ...)

A grade in the direct product of independent gradings, with one component per [`grading`](@ref). It is the product of
grades with different declared gradings, so that e.g. fermions graded by `Charge` and kets graded by `Bicyclic` can be
multiplied without defining `mul_grade(::Charge, ::Bicyclic)`, once both types declare their `grading`.

Products involving a `ProductGrade` are computed componentwise:
- components of the same grading are multiplied with [`mul_grade`](@ref), in order, so the product is associative;
- components that become `nothing` are dropped, and a single remaining component is returned unwrapped, so the
  product of the grades `Charge(1)`, `KET` and `Charge(-1)` is `KET`.

Components are matched by grading regardless of their order, so `f * k` and `k * f` have the same grade.
[`adjoint_grade`](@ref) acts componentwise. `mul_grade` and `adjoint_grade` within one grading must return a grade of
that grading or `nothing`, and this is checked. Construct a `ProductGrade` with at least two components, in any order,
with distinct declared gradings, none of them `nothing` or a `ProductGrade`.

# Example
```julia
# Fermions with spin, graded by charge and by spin, using `Charge` from the example in the README
struct Spin
    s::Int
end
spin(s) = iszero(s) ? nothing : Spin(s)
NonCommutativeProducts.grading(::Type{Spin}) = Spin # and likewise for Charge
NonCommutativeProducts.mul_grade(a::Spin, b::Spin) = spin(a.s + b.s)
NonCommutativeProducts.adjoint_grade(a::Spin) = spin(-a.s)
NonCommutativeProducts.grade(f::SpinFermion) = ProductGrade(Charge(f.dagger ? 1 : -1), Spin(f.up == f.dagger ? 1 : -1))
```
"""
struct ProductGrade
    grades::Tuple
    # Unchecked: the components must be valid, with distinct keys (see ProductGrade(gs...))
    ProductGrade(grades::Tuple, ::Val{:unchecked}) = new(grades)
end
function ProductGrade(gs...)
    length(gs) >= 2 || throw(ArgumentError("a ProductGrade needs at least two components, got $(length(gs))"))
    for g in gs
        (isnothing(g) || g isa ProductGrade) && throw(ArgumentError("invalid ProductGrade component $(repr(g))"))
        isnothing(grading(typeof(g))) && throw(ArgumentError("the grading of the ProductGrade component $(repr(g)) is \
            undeclared; define `NonCommutativeProducts.grading(::Type{$(typeof(g))})`"))
    end
    allunique(map(g -> grading(typeof(g)), gs)) || throw(ArgumentError("the components of a ProductGrade must have \
        different gradings, got $(join(map(repr, gs), ", "))"))
    return ProductGrade(gs, Val(:unchecked))
end

_grade_components(g::ProductGrade) = g.grades
_grade_components(g) = (g,)

# The component of `cs` with the key `key`, or nothing
function _component(cs::Tuple, key)
    for c in cs
        isequal(grading(typeof(c)), key) && return c
    end
    return nothing
end

# Components are compared by key, independently of their order
function Base.:(==)(g::ProductGrade, h::ProductGrade)
    length(g.grades) == length(h.grades) || return false
    for a in g.grades
        b = _component(h.grades, grading(typeof(a)))
        (isnothing(b) || !(a == b)) && return false
    end
    return true
end
function Base.isequal(g::ProductGrade, h::ProductGrade)
    length(g.grades) == length(h.grades) || return false
    for a in g.grades
        b = _component(h.grades, grading(typeof(a)))
        (isnothing(b) || !isequal(a, b)) && return false
    end
    return true
end
function Base.hash(g::ProductGrade, h::UInt)
    s = zero(UInt)
    for c in g.grades
        s += hash(c)
    end
    return hash(s, hash(ProductGrade, h))
end
# Sorted by key, so that the output is deterministic
function Base.show(io::IO, g::ProductGrade)
    cs = sort!(collect(Any, g.grades); by=c -> string(grading(typeof(c))))
    print(io, "ProductGrade(")
    join(io, (sprint(show, x; context=io) for x in cs), ", ")
    print(io, ")")
end

# A ProductGrade of nonidentity components with distinct keys, unwrapped if there are fewer than two
function _product_grade(cs::Vector{Any})
    isempty(cs) && return nothing
    length(cs) == 1 && return only(cs)
    return ProductGrade(Tuple(cs), Val(:unchecked))
end

# Merges the components of g and h by key, multiplying those of the same key with g's on the left
function _mul_product_grade(g, h)
    for x in (g, h)
        x isa ProductGrade || !isnothing(grading(typeof(x))) || throw(_undeclared_grading_error(g, h))
    end
    gs, hs = _grade_components(g), _grade_components(h)
    cs = Any[]
    for a in gs
        key = grading(typeof(a))
        b = _component(hs, key)
        if isnothing(b)
            push!(cs, a)
        else
            c = mul_grade(a, b)
            _is_closed(c, key) || throw(_closure_error(mul_grade, (a, b), c))
            isnothing(c) || push!(cs, c)
        end
    end
    for b in hs
        isnothing(_component(gs, grading(typeof(b)))) && push!(cs, b)
    end
    return _product_grade(cs)
end

# Componentwise
function _adjoint_grade(g::ProductGrade)
    cs = Any[]
    for a in g.grades
        c = adjoint_grade(a)
        _is_closed(c, grading(typeof(a))) || throw(_closure_error(adjoint_grade, (a,), c))
        isnothing(c) || push!(cs, c)
    end
    return _product_grade(cs)
end

# Inferred from the return type of grade, so this is a compile time constant for types without a grade method
_is_ungraded(::Type{S}) where {S} = Base.promote_op(grade, S) === Nothing
_keys_ungraded(::Type{<:NCMul{C,S}}) where {C,S} = _is_ungraded(S)
_keys_ungraded(::Type) = false

"""
    grades(x)

The set of grades of `x`. Grades never affect `==`, `hash` or any other property of the value.

For a nonzero `x`, these are the [`grade`](@ref)s of its terms, plus `nothing` if it has a nonzero scalar part.
A zero has the grades of how it was computed, from the grades `Γ(a)`, `Γ(b)` of the operands:

| operation | grades of a zero result |
|:--|:--|
| `a * b`, `a^n` | `Γ(a) ⊗ Γ(b) = {g ⋅ h : g ∈ Γ(a), h ∈ Γ(b)}` |
| `a + b`, `a - b` | `Γ(a) ∪ Γ(b)` |
| `a'` | `{adjoint_grade(g) : g ∈ Γ(a)}` |
| `α * a`, `zero(a)`, `zerovector(a)`, `copy`, `convert`, `sort` | `Γ(a)` |
| `mul!!(c, a, b)` | `Γ(c) ∪ Γ(a) ⊗ Γ(b)` |
| `add!!(a, b, α, β)`, `VectorInterface.add(!)` | `Γ(a) ∪ Γ(b)`, or `Γ(b)` if `β isa VectorInterface.Zero` |

Here `g ⋅ h` is the product of grades, [`mul_grade`](@ref) within one grading (see [`grading`](@ref)).

A zero without such a history (the number `0`, `zero(typeof(a))`, or a zero sum built from a `Dict` or by `ncmap`)
has no grades. The empty set is the unit of both `⊗` and `∪`.

A product without factors behaves like its prefactor, a number, so a zero one has no grades: `grades(2 * one(f))` is
`Set([nothing])`, but `grades(0 * one(f))` is empty. The table applies to products with factors and to sums.

# Example
```julia
f = Fermion(1) # with the charge grading in the example of `grade`
grades(f' * f + 1) == Set([nothing])
grades(f * f) == Set([-2]) # f * f == 0
grades(0 * (f' + f)) == Set([1, -1])
```
"""
grades(x) = _grades(x)
grades(x::NCAdd) = iszero(x) ? copy(_grades(x)) : _grades(x)

# Like grades, but may return the stored set of a zero sum, which must not be mutated
_grades(x) = Set{Any}((grade(x),))
_grades(x::Number) = iszero(x) ? Set{Any}() : Set{Any}((nothing,))
_grades(x::UniformScaling) = _grades(x.λ)
_grades(x::NCMul) = isempty(x.factors) ? _grades(prefactor(x)) : Set{Any}((grade(x),))
function _grades(x::NCAdd{C,K}) where {C,K}
    if iszero(x)
        gs = x.zero_grades
        return isnothing(gs) ? Set{Any}() : gs
    end
    _keys_ungraded(K) && return Set{Any}((nothing,))
    gs = Set{Any}(grade(k) for k in keys(x.dict))
    iszero(additive_coeff(x)) || push!(gs, nothing)
    return gs
end

"""
    grade_part(x, g)

The part of `x` with grade `g`, where the scalar part has grade `nothing`. `x` is the sum of `grade_part(x, g)` over
`g` in [`grades`](@ref)`(x)`.

# Example
```julia
f = Fermion(1) # with the charge grading in the example of `grade`
grade_part(f' * f + f + 1, nothing) == f' * f + 1
grade_part(f' * f + f + 1, -1) == f
```
"""
grade_part(x::Union{Number,UniformScaling}, g) = isnothing(g) ? x : zero(x)
# _grades rather than grade, which throws for a zero product without factors
grade_part(x, g) = g in _grades(x) ? x : grade_part(NCAdd(x), g)
function grade_part(x::NCAdd, g)
    coeff = isnothing(g) ? additive_coeff(x) : zero(additive_coeff(x))
    result = NCAdd(coeff, filter(kv -> isequal(grade(first(kv)), g), x.dict))
    return _with_zero_grades(result, () -> g in _grades(x) ? Set{Any}((g,)) : Set{Any}())
end

# `nothing` stands for ∅, as in `zero_grades`
_union_grades(::Nothing, gs) = gs
_union_grades(gs1, gs2) = union(gs1, gs2)
function _mul_grades(gs1, gs2)
    isempty(gs1) && return gs2
    isempty(gs2) && return gs1
    return Set{Any}(_mul_grade(g, h) for g in gs1 for h in gs2)
end
_adjoint_grades(gs) = Set{Any}(_adjoint_grade(g) for g in gs)

# The stored grades of `x` if it is zero, otherwise nothing (see the top of this file)
_zero_grades(x::NCAdd) = iszero(x) ? x.zero_grades : nothing

# Normalizes ∅ to `nothing`
_as_grades(::Nothing) = nothing
_as_grades(gs::Set{Any}) = isempty(gs) ? nothing : gs
_as_grades(gs) = isempty(gs) ? nothing : Set{Any}(gs)

function _set_zero_grades!(x::NCAdd, gs)
    x.zero_grades = _as_grades(gs)
    return x
end

# Restores the invariant at the exit of a public operation: stores `grades()` if `result` is zero, and clears stale
# grades otherwise. `grades` is only called for a zero result, behind an inference barrier so that compiled code
# doesn't depend on downstream grade methods (see _mul_effect in sorting.jl).
@inline function _with_zero_grades(result, grades)
    result isa NCAdd || return result
    if iszero(result)
        _set_zero_grades!(result, Base.inferencebarrier(grades)())
    elseif result.zero_grades !== nothing
        result.zero_grades = nothing
    end
    return result
end

"""
    check_grading(a, b)

Check that every term of [`mul_effect`](@ref)`(a, b)` has the product of the grades of `a` and `b` (see
[`grade`](@ref)). Zero terms are always allowed. Meant to be used in tests.

# Example
```julia
f = Fermion(1)
check_grading(f, f') # true
```
"""
function check_grading(a, b)
    # without autosort, as when sorting, so that the effect is checked as the sorter sees it
    effect = with(_autosort => false) do
        mul_effect(a, b)
    end
    isnothing(effect) && return true
    return _has_grade(effect, _mul_grade(grade(a), grade(b)), _mul_grade(grade(b), grade(a)))
end
# `swapgrade` is the grade of b * a, the product that Swap stands for
_has_grade(x::Number, g, swapgrade) = iszero(x) || isnothing(g)
_has_grade(x::NCMul, g, swapgrade) = iszero(x) || isequal(grade(x), g)
_has_grade(x::NCAdd, g, swapgrade) = iszero(x) || all(isequal(g), _grades(x))
_has_grade(x::Swap, g, swapgrade) = iszero(x.λ) || isequal(swapgrade, g)
_has_grade(x::AddTerms, g, swapgrade) = all(term -> _has_grade(term, g, swapgrade), x.terms)
_has_grade(x, g, swapgrade) = isequal(grade(x), g)

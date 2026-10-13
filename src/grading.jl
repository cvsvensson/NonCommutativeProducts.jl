# Gradings. A grading assigns a grade to every atom, and the grade of a product is the product of the grades of its
# factors. `nothing` is the identity grade, which nonzero numbers and ungraded atoms have. Below, Γ(x) = grades(x) is
# a set of grades, and ∅ is the empty set. The contract:
#
# 1. Value. ==, isequal, hash, iszero, isscalar and scalar depend only on the value, never on the representation or
#    on the grades. All zeros are equal, to each other and to 0.
# 2. The grades of a nonzero expression are exact: the grades of its homogeneous parts, with `nothing` for a nonzero
#    scalar part. They are computed from the terms.
# 3. The grades of a zero come from how it was computed, applied once per operation to the grades of the operands:
#    Γ(a)⊗Γ(b) = {mul_grade(g, h)} for a product (also ^ and sorting a product), Γ(a) ∪ Γ(b) for a sum, the adjoint
#    grades for an adjoint, Γ(x) for zero, zerovector(!), copy, convert, anyadd and sort, Γ(α)⊗Γ(x) = Γ(x) for
#    scaling, Γ(c) ∪ Γ(a)⊗Γ(b) for mul!!(c, a, b), and Γ(a) ∪ Γ(b) for add!!(a, b, α, β), or Γ(b) if β isa Zero.
#    The destination of an in-place operation counts as an operand only if its old value enters the result, so what a
#    mutated object held before never affects the grades of the zero it now holds.
# 4. A bare zero, i.e. one without a derivation (a zero number, a zero product without factors, a zero sum of unknown
#    origin), has Γ = ∅, and ∅ is the unit of both ⊗ and ∪: a bare zero takes on the grades of what it is combined
#    with. This is never false, since zero lies in every graded piece.
# 5. grade(x) == only(grades(x)), and it throws an ArgumentError if there are no grades or several. The only primitive
#    is the user's grade(atom).
#
# A zero sum has no terms, so its grades are stored in the `zero_grades` field of NCAdd, where `nothing` stands for ∅
# (an empty set is never stored). Rule 3 is enforced through invariant Z on every NCAdd that a public operation
# returns or mutates: if it is zero, `zero_grades` holds the grades of rule 3 for the last public operation that wrote
# it, and if it is nonzero, `zero_grades === nothing`. Public operations restore Z at their exit with
# _with_zero_grades, which writes the grades of a zero result and clears them from a nonzero one. Internal helpers
# (set_coeff!, add_to_coeff!, filter_ncadd!!, _set_ncadd!, _add_scaled!!, bubble_sort!, ...) never maintain the
# grades, and may leave them stale for their caller to restore. Editing the fields directly, or the constructor flags
# filter_zeros=false and filter_scalars=false, may leave the grades unknown (∅), but never wrong.
#
# The rules agree with the values only if every mul_effect respects the grading, which check_grading checks. Under that
# assumption, add!! and mul!! may skip the grades of a nonzero destination when the result is zero: if c + v == 0 with
# c != 0, then c == -v, so the grades of c are among those of v.

"""
    grade(atom)

The grade of `atom` in a grading of the algebra. Define it for your types to make them graded. The default is
`nothing`, the identity grade, so that types without a method are ungraded. Grades are multiplied with
[`mul_grade`](@ref), and every term that [`mul_effect`](@ref)`(a, b)` produces must have the grade
`mul_grade(grade(a), grade(b))`, which [`check_grading`](@ref) checks.

For any other expression, `grade` is derived: it is the only element of [`grades`](@ref)`(x)`, and it throws an
`ArgumentError` if there are several grades, or none. So a nonzero number has the grade `nothing`, and a product with
factors has the product of the grades of its factors, regardless of the coefficient (`grade(0 * a) == grade(a)`).
A zero without grades, such as the number `0` or a zero sum of unknown origin, has no grade, and `grade` throws.

# Example
```julia
# Fermions with charge as the grade, so that number-conserving operators have the identity grade
NonCommutativeProducts.grade(f::Fermion) = f.creation ? 1 : -1
NonCommutativeProducts.mul_grade(g::Int, h::Int) = iszero(g + h) ? nothing : g + h
NonCommutativeProducts.adjoint_grade(g::Int) = -g
```
"""
grade(x) = nothing
grade(x::Number) = iszero(x) ? throw(_no_grade_error(x, Set{Any}())) : nothing
grade(x::UniformScaling) = grade(x.λ)
function grade(x::NCMul{C,S}) where {C,S}
    # a product without factors is a number, and if it is zero it has no grades
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

The product of the grades `g` and `h`, which is only called when neither of them is `nothing`. Define it for the
grades returned by [`grade`](@ref). It must be associative, and return `nothing` for the identity grade.
"""
function mul_grade end

"""
    adjoint_grade(g)

The grade of the adjoint of an expression of grade `g`, which must agree with `grade(a')` for every atom `a` of grade
`g`. It is only called for grades other than `nothing`, and only needed for adjoints of zero sums.
"""
function adjoint_grade end

_mul_grade(::Nothing, ::Nothing) = nothing
_mul_grade(::Nothing, h) = h
_mul_grade(g, ::Nothing) = g
_mul_grade(g, h) = mul_grade(g, h)
_adjoint_grade(::Nothing) = nothing
_adjoint_grade(g) = adjoint_grade(g)

# Whether every atom of type S has the identity grade. This is decided from the inferred return type of grade, so
# for a concrete type without a grade method it is a compile time constant, and no grades are computed.
_is_ungraded(::Type{S}) where {S} = Base.promote_op(grade, S) === Nothing
_keys_ungraded(::Type{<:NCMul{C,S}}) where {C,S} = _is_ungraded(S)
_keys_ungraded(::Type) = false

"""
    grades(x)

The set of grades of `x`. Grades are metadata: they never affect `==`, `hash` or any other property of the value.

For a nonzero `x` they are exact: the [`grade`](@ref) of every term, and `nothing` if `x` has a nonzero scalar part.
So a product or an atom has the set with its grade, and a nonzero number has the identity grade `nothing`.

All zeros are equal, but a zero still has the grades of how it was computed. A zero product with factors has its
grade (e.g. `grades(0 * a) == grades(a)`), and a zero sum has the grades given by the operation that produced it, from
the grades `Γ(a)`, `Γ(b)` of its operands:

| operation | grades of a zero result |
|:--|:--|
| `a * b`, `a^n` | `Γ(a) ⊗ Γ(b) = {mul_grade(g, h) : g ∈ Γ(a), h ∈ Γ(b)}` |
| `a + b`, `a - b` | `Γ(a) ∪ Γ(b)` |
| `a'` | `{adjoint_grade(g) : g ∈ Γ(a)}` |
| `α * a`, `zero(a)`, `zerovector(a)`, `copy`, `convert`, `sort` | `Γ(a)` |
| `mul!!(c, a, b)` | `Γ(c) ∪ Γ(a) ⊗ Γ(b)` |
| `add!!(a, b, α, β)`, `VectorInterface.add(!)` | `Γ(a) ∪ Γ(b)`, or `Γ(b)` if `β isa VectorInterface.Zero` |

An in-place operation follows the same rule, so what the destination held before only matters if its old value
enters the result.

A bare zero, i.e. one without a derivation, has no grades (∅): the number `0`, a zero product without factors such as
`zero(typeof(a))`, and a zero sum of unknown origin, e.g. one built from a `Dict` or by `ncmap`. ∅ is the unit of both
`⊗` and `∪`, so a bare zero takes on the grades of whatever it is combined with, e.g. `grades(0 + (b - b)) == grades(b)`.

# Example
```julia
f = Fermion(1) # with the charge grading in the example of `grade`, and autosort enabled
grades(f' * f + 1) == Set([nothing])
grades(f * f) == Set([-2]) # f * f == 0
grades(0 * (f' + f)) == Set([1, -1])
```
"""
grades(x) = _grades(x)
grades(x::NCAdd) = iszero(x) ? copy(_grades(x)) : _grades(x)

# Like grades, but the result may be shared with the zero sum it is read from, so it must not be mutated.
# Grade sets are always a Set{Any}, and are never mutated once they are stored in an NCAdd.
_grades(x) = Set{Any}((grade(x),))
_grades(x::Number) = iszero(x) ? Set{Any}() : Set{Any}((nothing,))
_grades(x::UniformScaling) = _grades(x.λ)
# a product without factors is a number, so if it is zero it is a bare zero
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

The part of `x` with grade `g`. For a sum, it is the sum of the terms of grade `g`, where the scalar part has the
grade `nothing`, so that `x` is the sum of `grade_part(x, g)` over `g` in [`grades`](@ref)`(x)`. A product, an atom
or a number is returned as it is if its [`grade`](@ref) is `g`, and is otherwise replaced by a zero. A zero result
has the grades `grades(x) ∩ {g}`.

# Example
```julia
f = Fermion(1) # with the charge grading in the example of `grade`, and autosort enabled
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

# The operations on grade sets in rule 3. ∅ is the unit of both ∪ and ⊗ (rule 4), and `nothing` stands for ∅, as in
# the `zero_grades` field, so an operand whose grades are `nothing` is skipped. The result may be one of the arguments,
# which is fine since grade sets are never mutated.
_union_grades(::Nothing, gs) = gs
_union_grades(gs1, gs2) = union(gs1, gs2)
function _mul_grades(gs1, gs2)
    isempty(gs1) && return gs2
    isempty(gs2) && return gs1
    return Set{Any}(_mul_grade(g, h) for g in gs1 for h in gs2)
end
_adjoint_grades(gs) = Set{Any}(_adjoint_grade(g) for g in gs)

# The stored grades of `x` if it is zero, otherwise nothing. Operations that mutate `x` read this first. The grades of
# a nonzero destination are skipped, which is sound under the assumption that mul_effect respects the grading (see
# the top of this file).
_zero_grades(x::NCAdd) = iszero(x) ? x.zero_grades : nothing

# Grades as stored in `zero_grades`: a nonempty Set{Any}, or `nothing` for ∅, so that ∅ has a single representation
_as_grades(::Nothing) = nothing
_as_grades(gs::Set{Any}) = isempty(gs) ? nothing : gs
_as_grades(gs) = isempty(gs) ? nothing : Set{Any}(gs)

function _set_zero_grades!(x::NCAdd, gs)
    x.zero_grades = _as_grades(gs)
    return x
end

# Restores invariant Z at the exit of a public operation (see the top of this file). If `result` is a zero sum, store
# `grades()` as its grades, and if it is a nonzero sum, clear any grades it still holds from an earlier value.
# `grades` is a function of the operands, which must still be intact, and is only called when the result is zero, so
# the nonzero path costs one iszero check and a field read. It is called behind an inference barrier, so that the
# compiled code of the operations doesn't depend on the grade methods defined downstream (see _mul_effect in
# sorting.jl). In nested public operations the outermost call writes last, so the inner calls are harmless.
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

Check that every term of [`mul_effect`](@ref)`(a, b)` has the grade `mul_grade(grade(a), grade(b))`, as the grading
requires. The terms of `Swap(λ)` (the product `λ * b * a`) and of `AddTerms` are checked individually. Zero terms and
an effect of `nothing` are always allowed. Returns a `Bool`, and is meant to be called in the tests of the rules.

# Example
```julia
f = Fermion(1)
check_grading(f, f') # true
```
"""
function check_grading(a, b)
    effect = mul_effect(a, b)
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

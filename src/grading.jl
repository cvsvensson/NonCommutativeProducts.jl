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
Every term of [`mul_effect`](@ref)`(a, b)` must have the grade [`mul_grade`](@ref)`(grade(a), grade(b))`, see
[`check_grading`](@ref).

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

The product of the grades `g` and `h`, neither of which is `nothing`. It must be associative, and return `nothing` for
the identity grade.
"""
function mul_grade end

"""
    adjoint_grade(g)

The grade of the adjoint of an expression of grade `g` (not `nothing`). Only needed for adjoints of zero sums.
"""
function adjoint_grade end

_mul_grade(::Nothing, ::Nothing) = nothing
_mul_grade(::Nothing, h) = h
_mul_grade(g, ::Nothing) = g
_mul_grade(g, h) = mul_grade(g, h)
_adjoint_grade(::Nothing) = nothing
_adjoint_grade(g) = adjoint_grade(g)

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
| `a * b`, `a^n` | `Γ(a) ⊗ Γ(b) = {mul_grade(g, h) : g ∈ Γ(a), h ∈ Γ(b)}` |
| `a + b`, `a - b` | `Γ(a) ∪ Γ(b)` |
| `a'` | `{adjoint_grade(g) : g ∈ Γ(a)}` |
| `α * a`, `zero(a)`, `zerovector(a)`, `copy`, `convert`, `sort` | `Γ(a)` |
| `mul!!(c, a, b)` | `Γ(c) ∪ Γ(a) ⊗ Γ(b)` |
| `add!!(a, b, α, β)`, `VectorInterface.add(!)` | `Γ(a) ∪ Γ(b)`, or `Γ(b)` if `β isa VectorInterface.Zero` |

A zero without such a history (the number `0`, `zero(typeof(a))`, or a zero sum built from a `Dict` or by `ncmap`)
has no grades. The empty set is the unit of both `⊗` and `∪`.

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

Check that every term of [`mul_effect`](@ref)`(a, b)` has the grade `mul_grade(grade(a), grade(b))`. Zero terms are
always allowed. Meant to be used in tests.

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

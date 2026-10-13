# Gradings. A grading assigns a grade to every atom, and the grade of a product is the product of the grades of its
# factors. `nothing` is the identity grade, which numbers and ungraded atoms have.
#
# The grades of a nonzero expression are derived from its terms. A zero sum has no terms, so its grades are instead
# stored in the `zero_grades` field of NCAdd when an operation produces it, from the grades of the operands (see
# _with_zero_grades). 

"""
    grade(atom)

The grade of `atom` in a grading of the algebra. Define it for your types to make them graded. The default is
`nothing`, the identity grade, so that types without a method are ungraded. Grades are multiplied with
[`mul_grade`](@ref), and every term that [`mul_effect`](@ref)`(a, b)` produces must have the grade
`mul_grade(grade(a), grade(b))`, which [`check_grading`](@ref) checks.

For a product, `grade` is the product of the grades of its factors, regardless of the coefficient. For a sum, it is
the only element of [`grades`](@ref), and it throws if the sum has several grades or none.

# Example
```julia
# Fermions with charge as the grade, so that number-conserving operators have the identity grade
NonCommutativeProducts.grade(f::Fermion) = f.creation ? 1 : -1
NonCommutativeProducts.mul_grade(g::Int, h::Int) = iszero(g + h) ? nothing : g + h
NonCommutativeProducts.adjoint_grade(g::Int) = -g
```
"""
grade(x) = nothing
grade(::Number) = nothing
grade(::UniformScaling) = nothing
function grade(x::NCMul{C,S}) where {C,S}
    _is_ungraded(S) && return nothing
    return foldl((g, f) -> _mul_grade(g, grade(f)), x.factors; init=nothing)
end
grade(x::NCAdd) = only(grades(x))

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

The set of grades of the terms of `x`. For a product or an atom, it is the set with its [`grade`](@ref). For a nonzero
sum, it contains the grade of every term, and `nothing` if the sum has a nonzero scalar part. A nonzero number has
the identity grade `nothing`, and `0` has no grades.

A zero sum has no terms, so its grades are instead those of the operands it was computed from: the products of their
grades for a product, the union for a sum, and the grades of the operand for scaling, `zero(x)`, `zerovector`,
`copy`, `convert` and sorting. A zero sum whose origin is unknown, e.g. one built from a `Dict` or by `ncmap`, has no
grades.

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

_union_grades(::Nothing, gs) = gs
_union_grades(gs1, gs2) = union(gs1, gs2)
_mul_grades(gs1, gs2) = Set{Any}(_mul_grade(g, h) for g in gs1 for h in gs2)
_adjoint_grades(gs) = Set{Any}(_adjoint_grade(g) for g in gs)

# The stored grades of `x` if it is zero, otherwise nothing. Operations that mutate `x` read this first.
_zero_grades(x::NCAdd) = iszero(x) ? x.zero_grades : nothing

_as_grades(gs::Union{Nothing,Set{Any}}) = gs
_as_grades(gs) = Set{Any}(gs)

function _set_zero_grades!(x::NCAdd, gs)
    x.zero_grades = _as_grades(gs)
    return x
end

# The typing rule for zero results: if `result` is a zero sum, store `grades()` as its grades. `grades` is a function
# of the operands, which must still be intact, and is only called when the result is zero, so the nonzero path costs
# one iszero check. It is called behind an inference barrier, so that the compiled code of the operations doesn't
# depend on the grade methods defined downstream (see _mul_effect in sorting.jl).
@inline function _with_zero_grades(result, grades)
    result isa NCAdd && iszero(result) || return result
    return _set_zero_grades!(result, Base.inferencebarrier(grades)())
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

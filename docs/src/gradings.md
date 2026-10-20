# Gradings

```@meta
DocTestSetup = quote
    using KrylovKit, LinearAlgebra
    using NonCommutativeProducts
    import NonCommutativeProducts as NC
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
    struct Ket
        n::Int
        bra::Bool
        Ket(n::Int, bra::Bool=false) = new(n, bra)
    end
    Base.show(io::IO, k::Ket) = print(io, k.bra ? "⟨" : "|", k.n, k.bra ? "|" : "⟩")
    Base.adjoint(k::Ket) = Ket(k.n, !k.bra)
    NC.@nc Ket
    NC.mul_effect(a::Ket, b::Ket) = a.bra && !b.bra ? Int(a.n == b.n) : nothing
    k0, k1 = Ket(0), Ket(1)
    NC.enable_autosort!()
end
```

The examples on this page use the `Fermion` type from the [Home](@ref "Example: Fermions") page and the `Ket` type from the [Linear algebra](@ref) page, with autosort enabled.

A grading gives every factor a grade, so that the grade of a product is the product of the grades of its factors. Gradings can be used to keep track of what kind of object an expression is, e.g. a ket, a bra or an operator, or of quantities such as particle number. The grade `nothing` is the identity grade, which numbers and types without a grade have. Define these for your types:
* [`grade(a)`](@ref NonCommutativeProducts.grade): the grade of the factor `a`.
* [`mul_grade(g, h)`](@ref NonCommutativeProducts.mul_grade): the product of two grades of one grading. It must be associative, return a grade of the same grading, and return `nothing` for the identity.
* [`adjoint_grade(g)`](@ref NonCommutativeProducts.adjoint_grade): the grade of `a'` when `a` has grade `g`.

Then [`grades(x)`](@ref NonCommutativeProducts.grades) returns the set of grades of the terms of `x`, and `grade(x)` returns the single grade of `x`, throwing an error if there isn't exactly one.

The value `0` doesn't say what kind of zero it is. Instead, a zero keeps the grades of the expressions it was derived from. For the kets of the [Linear algebra](@ref) page, a ket has grade (1,0) and a bra (0,1).
```jldoctest gradings
julia> struct Bicyclic
           m::Int
           n::Int
       end

julia> NC.grade(k::Ket) = k.bra ? Bicyclic(0, 1) : Bicyclic(1, 0)

julia> function NC.mul_grade(g::Bicyclic, h::Bicyclic)
           c = max(g.n, h.m)
           r = Bicyclic(g.m - g.n + c, h.n - h.m + c)
           return r == Bicyclic(0, 0) ? nothing : r
       end

julia> NC.adjoint_grade(g::Bicyclic) = Bicyclic(g.n, g.m)

julia> z = k0 * k1' * k0 # an annihilated ket
0I

julia> z == 0
true

julia> NC.grades(z) # a zero ket
Set{Any} with 1 element:
  Bicyclic(1, 0)

julia> NC.grades(z') # a zero bra
Set{Any} with 1 element:
  Bicyclic(0, 1)
```

As another example, give the fermions a charge, to check that a Hamiltonian conserves the number of particles:
```jldoctest gradings
julia> struct Charge
           q::Int
       end

julia> NC.grade(f::Fermion) = Charge(f.dagger ? 1 : -1)

julia> NC.mul_grade(a::Charge, b::Charge) = iszero(a.q + b.q) ? nothing : Charge(a.q + b.q)

julia> NC.adjoint_grade(a::Charge) = Charge(-a.q)

julia> H = sum(Fermion(n)' * Fermion(n + 1) + Fermion(n + 1)' * Fermion(n) for n in 1:3) + Fermion(1)' * Fermion(1);

julia> NC.grades(H) # H conserves the particle number
Set{Any} with 1 element:
  nothing

julia> NC.grades(Fermion(1)' + Fermion(2)' * Fermion(1)' * Fermion(1))
Set{Any} with 1 element:
  Charge(1)
```
To mix grades of different types, declare which grading each type belongs to with [`grading`](@ref NonCommutativeProducts.grading):
```jldoctest gradings
julia> NC.grading(::Type{Charge}) = Charge

julia> NC.grading(::Type{Bicyclic}) = Bicyclic

julia> NC.@commutative Fermion Ket
```
Grades of the same type, or of different types with the same declared grading, are multiplied with `mul_grade`. Grades of different declared gradings are independent and combine into a [`NC.ProductGrade`](@ref NonCommutativeProducts.ProductGrade), with one component per grading. Mixing grades where either type has no declared grading throws an `ArgumentError`. So `Fermion(1)' * k0` has the grade `ProductGrade(Bicyclic(1, 0), Charge(1))` in either order.
```jldoctest gradings
julia> NC.grade(Fermion(1)' * k0)
ProductGrade(Bicyclic(1, 0), Charge(1))

julia> NC.grade(k0 * Fermion(1)')
ProductGrade(Bicyclic(1, 0), Charge(1))
```
Components of the same grading are multiplied with each other, and components that become `nothing` are dropped, so `Fermion(1)' * k0 * Fermion(1)` has the grade `Bicyclic(1, 0)`.
```jldoctest gradings
julia> NC.grade(Fermion(1)' * k0 * Fermion(1))
Bicyclic(1, 0)
```
Several types can share a grading, e.g. `NC.grading(::Type{<:U1}) = U1`, and then you define `mul_grade` between them. A type can also carry several gradings, e.g. `NC.grade(f::SpinfulFermion) = NC.ProductGrade(Charge(1), Spin(1))`, if `Charge` and `Spin` have declared gradings.

Every term that `mul_effect(a, b)` returns must have the product of the grades of `a` and `b`. Use [`check_grading(a, b)`](@ref NonCommutativeProducts.check_grading) in your tests to check this:
```jldoctest gradings
julia> NC.check_grading(Fermion(1), Fermion(1)')
true
```

```@meta
DocTestSetup = nothing
```

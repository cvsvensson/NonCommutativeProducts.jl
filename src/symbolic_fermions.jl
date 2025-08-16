struct NormalOrdering <: AbstractOrdering end


struct SymbolicFermionBasis
    name::Symbol
    universe::UInt64
end
Base.hash(x::SymbolicFermionBasis, h::UInt) = hash(x.name, hash(x.universe, h))

"""
    @fermions a b ...

Create one or more fermion species with the given names. Indexing into fermions species
gives a concrete fermion. Fermions in one `@fermions` block anticommute with each other, 
and commute with fermions in other `@fermions` blocks.

# Examples:
- `@fermions a b` creates two species of fermions that anticommute:
    - `a[1]' * a[1] + a[1] * a[1]' == 1`
    - `a[1]' * b[1] + b[1] * a[1]' == 0`
- `@fermions a; @fermions b` creates two species of fermions that commute with each other:
    - `a[1]' * a[1] + a[1] * a[1]' == 1`
    - `a[1] * b[1] - b[1] * a[1] == 0`

See also [`@majoranas`](@ref)
"""
macro fermions(xs...)
    universe = hash(xs)
    defs = map(xs) do x
        :($(esc(x)) = SymbolicFermionBasis($(Expr(:quote, x)), $universe))
    end
    Expr(:block, defs...,
        :(tuple($(map(x -> esc(x), xs)...))))
end
Base.:(==)(a::SymbolicFermionBasis, b::SymbolicFermionBasis) = a.name == b.name && a.universe == b.universe
Base.getindex(f::SymbolicFermionBasis, is...) = Fermion(false, is, f)
Base.getindex(f::SymbolicFermionBasis, i) = Fermion(false, i, f)

struct Fermion{L,B}
    creation::Bool
    label::L
    basis::B
end
Base.adjoint(x::Fermion) = Fermion(!x.creation, x.label, x.basis)
Base.iszero(x::Fermion) = false
function Base.show(io::IO, x::Fermion)
    print(io, x.basis.name, x.creation ? "†" : "")
    if Base.isiterable(typeof(x.label))
        Base.show_delim_array(io, x.label, "[", ",", "]", false)
    else
        print(io, "[", x.label, "]")
    end
end
function Base.isless(a::Fermion, b::Fermion)
    if a.basis.universe !== b.basis.universe
        a.basis.universe < b.basis.universe
    elseif a.creation == b.creation
        a.basis.name == b.basis.name && return a.label < b.label
        a.basis.name < b.basis.name
    else
        a.creation > b.creation
    end
end
Base.:(==)(a::Fermion, b::Fermion) = a.creation == b.creation && a.label == b.label && a.basis == b.basis
Base.hash(a::Fermion, h::UInt) = hash(a.creation, hash(a.label, hash(a.basis, h)))

# ordered_product(as::NCMul, bs::NCMul, ::NormalOrdering) = canonicalize!(normal_order(ordered_product(as, bs, NaiveOrdering())))

TermInterface.head(::T) where {T<:Fermion} = T
TermInterface.iscall(::Fermion) = true
TermInterface.isexpr(::Fermion) = true
TermInterface.maketerm(::Type{Q}, head::Type{T}, args, metadata) where {Q<:Union{Fermion,<:NCMul,<:NCAdd},T<:Fermion} = T(args...)


function ordered_product(a::Fermion, b::Fermion, ::NormalOrdering)
    a_uni = a.basis.universe
    b_uni = b.basis.universe
    if a == b
        0
    elseif a < b
        NCMul(1, [a, b])
    elseif a > b
        NCMul((-1)^(a_uni == b_uni), [b, a]) + Int(a.label == b.label && a.basis == b.basis)
    else
        throw(ArgumentError("Don't know how to multiply $a * $b"))
    end
end
function Base.:^(a::Fermion, b)
    if b isa Number && iszero(b)
        1
    elseif b isa Number && b == 1
        a
    elseif b isa Integer && b >= 2
        0
    else
        throw(ArgumentError("Invalid exponent $b"))
    end
end

@testitem "SymbolicFermions" begin
    using Symbolics, LinearAlgebra
    @fermions f c
    @fermions b
    @variables a::Real z::Complex
    f1 = f[:a]
    f2 = f[:b]
    f3 = f[1, :↑]

    # Test canonical commutation relations
    @test f1' * f1 + f1 * f1' == 1
    @test iszero(f1 * f2 + f2 * f1)
    @test iszero(f1' * f2 + f2 * f1')

    # c anticommutes with f
    @test iszero(f1' * c[1] + c[1] * f1')
    # b commutes with f
    @test iszero(f1' * b[1] - b[1] * f1')

    @test_nowarn display(f1)
    @test_nowarn display(f3)
    @test_nowarn display(1 * f1)
    @test_nowarn display(2 * f3)
    @test_nowarn display(1 + f1)
    @test_nowarn display(1 + f3)
    @test_nowarn display(1 + a * f2 - 5 * f1 + 2 * z * f1 * f2)

    @test iszero(f1 - f1)
    @test iszero(f1 * f1)
    @test f1 * f2 isa NonCommutativeProducts.NCMul
    @test iszero(2 * f1 - 2 * f1)
    @test iszero(0 * f1)
    @test 2 * f1 isa NonCommutativeProducts.NCMul
    @test iszero(f1 * 0)
    @test iszero(f1^2)
    @test iszero(0 * (f1 + f2))
    @test iszero((f1 + f2) * 0)
    @test iszero(f1 * f2 * f1)
    f12 = f1 * f2
    @test iszero(f12'' - f12)
    @test iszero(f12 * f12)
    @test iszero(f12' * f12')
    nf1 = f1' * f1
    @test nf1^2 == nf1
    @test f1' * f1 isa NonCommutativeProducts.NCMul
    @test f1 * f1' isa NonCommutativeProducts.NCAdd

    @test 1 + (f1 + f2) == 1 + f1 + f2 == f1 + f2 + 1 == f1 + 1 + f2 == 1 * f1 + f2 + 1 == f1 + 0.5 * f2 + 1 + (0 * f1 + 0.5 * f2) == (0.5 + 0.5 * f1 + 0.2 * f2) + 0.5 + (0.5 * f1 + 0.8 * f2) == (1 + f1' + (1 * f2)')'
    @test iszero((2 * f1) * (2 * f1))
    @test iszero((2 * f1)^2)
    @test (2 * f2) * (2 * f1) == -4 * f1 * f2
    @test f1 == (f1 * (f1 + 1)) == (f1 + 1) * f1
    @test iszero(f1 * (f1 + f2) * f1)
    @test (f1 * (f1 + f2)) == f1 * f2
    @test (2nf1 - 1) * (2nf1 - 1) == 1

    @test (1 * f1) * f2 == f1 * f2
    @test (1 * f1) * (1 * f2) == f1 * f2
    @test f1 * f2 == f1 * (1 * f2) == f1 * f2
    @test f1 - 1 == (1 * f1) - 1 == (0.5 + f1) - 1.5

    for op in [f1, f1 + f2, f1' * f2, f1 * f2 + 1]
        @test op + 1.0I == op + 1.0
        @test op - 1.0I == op - 1.0 == -(1.0I - op) == -(1.0 - op)
    end

    ex = 2 * f1
    @test NonCommutativeProducts.head(ex) == (*)
    @test NonCommutativeProducts.children(ex) == [2, ex.factors...]
    @test NonCommutativeProducts.operation(ex) == (*)
    @test NonCommutativeProducts.arguments(ex) == [2, ex.factors...]
    @test NonCommutativeProducts.isexpr(ex)
    @test NonCommutativeProducts.iscall(ex)
    ex = 2 * f1 + 1
    @test NonCommutativeProducts.head(ex) == (+)
    @test NonCommutativeProducts.children(ex) == [1, 2 * f1]
    @test NonCommutativeProducts.operation(ex) == (+)
    @test NonCommutativeProducts.arguments(ex) == [1, 2 * f1]
    @test NonCommutativeProducts.isexpr(ex)
    @test NonCommutativeProducts.iscall(ex)

    ex = f1
    @test NonCommutativeProducts.head(ex) <: NonCommutativeProducts.Fermion
    @test NonCommutativeProducts.children(ex) == [false, :a, f]
    @test NonCommutativeProducts.operation(ex) == NonCommutativeProducts.Fermion
    @test NonCommutativeProducts.arguments(ex) == [false, :a, f]
    @test NonCommutativeProducts.isexpr(ex)
    @test NonCommutativeProducts.iscall(ex)

    @test substitute(f1, f1 => f2) == f2
    @test substitute(f1', f1' => f2) == f2
    @test substitute(f1', f1 => f2) == f1'
    @test substitute(f1 + f2, f1 => f2) == 2 * f2

    @test substitute(2 * f1, 2 => 3) == 3 * f1
    @test iszero(substitute(a * f1 + 1, a => a^2) - (a^2 * f1 + 1))
    @test iszero(substitute(a * f1 * f2 - f1 + 1 + 0.5 * f2' * f2, f1 => f2) - (a * f2 * f2 - f2 + 1 + 0.5 * f2' * f2))
    @test iszero(substitute(a * f1 + a + a * f1 * f2 * f1', a => 0))

    @test substitute(f[1], 1 => 2) == f[2]
    @test substitute(f[:a]' * f[:b] + 1, :a => :b) == f[:b]' * f[:b] + 1
end

TermInterface.operation(::Fermion) = Fermion
TermInterface.arguments(a::Fermion) = [a.creation, a.label, a.basis]
TermInterface.children(a::Fermion) = arguments(a)
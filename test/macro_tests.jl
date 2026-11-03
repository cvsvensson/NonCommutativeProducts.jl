@testitem "@nc rejects what it can't register" begin
    import NonCommutativeProducts: @nc
    struct NN <: Number end
    struct NotNumber end
    abstract type AbstractAtom end
    struct NCAtom <: AbstractAtom end

    @test_throws ArgumentError @nc Int
    @test_throws ArgumentError @nc NN
    @test_throws ArgumentError @nc Real
    @test_throws ArgumentError @nc Union{NotNumber,Int}
    @test_throws ArgumentError @nc NotNumber, NCAtom
    # nothing is defined for the type that came before a rejected one
    @test_throws ArgumentError @nc NotNumber Int
    @test_throws MethodError NotNumber() + NotNumber()

    # abstract, parametric and Union types that hold no numbers, and repeated registration
    @nc AbstractAtom
    @nc NCAtom
    @test NCAtom() + NCAtom() == 2 * NCAtom()
    struct P{T} end
    struct U1 end
    struct U2 end
    @nc P
    @nc Union{U1,U2}
    @test P{Int}() + P{Int}() == 2 * P{Int}()
    @test U1() + U2() == U2() + U1()
end

@testitem "@commutative checks its types" begin
    import NonCommutativeProducts: @nc, @commutative
    struct C1 id::Int end
    struct C2 end
    struct C3 end
    struct Unregistered end
    abstract type AbstractC end
    struct C4 <: AbstractC end

    @test_throws ArgumentError @commutative C2 Unregistered
    @nc C1 C2 C3 AbstractC

    @test_throws ArgumentError @commutative C1 C1
    @test_throws ArgumentError @commutative AbstractC C4
    @test_throws ArgumentError @commutative C1, C2

    @commutative C1 C2
    @commutative C1 C2 # repeating the same order is fine
    @test C2() * C1(1) == C1(1) * C2()
    @test_throws ArgumentError @commutative C2 C1

    @commutative Union{C1,C2} C3
    @test C3() * C2() == C2() * C3()
end

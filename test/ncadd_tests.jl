@testitem "NCAdd: failed in-place updates leave the target unchanged" setup = [Fermions] begin
    using VectorInterface
    import NonCommutativeProducts as NC
    f1, f2 = Fermion(1), Fermion(2)

    # Int coefficients scaled by 0.5: the coefficient 2 * 0.5 fits in Int, the term f1 * 0.5 doesn't
    x = 2 + f1 + 2 * f2
    @test_throws InexactError VectorInterface.scale!(x, 0.5)
    @test x == 2 + f1 + 2 * f2
    @test_throws InexactError NC.scale!(x, 0.5)
    @test x == 2 + f1 + 2 * f2
    @test_throws InexactError NC.add!(x, 1, 1, 0.5)
    @test x == 2 + f1 + 2 * f2
    # coefficient that doesn't fit
    @test_throws InexactError NC.add!(x, 0.5)
    @test x == 2 + f1 + 2 * f2

    y = 2 + f1
    @test_throws InexactError VectorInterface.scale!(y, 2 + f1 + 2 * f2, 0.5)
    @test y == 2 + f1
    @test_throws InexactError VectorInterface.scale!(y, NC.NCMul(1, [f2]), 0.5)
    @test y == 2 + f1
    @test_throws InexactError VectorInterface.add!(y, 2 + f1 + 2 * f2, 0.5, 1)
    @test y == 2 + f1
    @test_throws InexactError VectorInterface.add!(y, 2 + 2 * f2, 1, 0.5)
    @test y == 2 + f1

    # and when they succeed, they mutate
    y = 2 + f1
    @test VectorInterface.scale!(y, 2 + f1 + 2 * f2, 2) === y
    @test y == 4 + 2 * f1 + 4 * f2
    @test VectorInterface.scale!(y, 3 * f2, 2) === y
    @test y == 6 * f2
    @test VectorInterface.add!(y, 2 + f1, 1, 2) === y
    @test y == 2 + f1 + 12 * f2
end

@testitem "NCAdd: set_coeff!!" setup = [Fermions] begin
    import NonCommutativeProducts as NC
    f1 = Fermion(1)
    a = f1 + 1
    b = NC.set_coeff!!(a, 3)
    @test b === a
    @test a == f1 + 3
    c = NC.set_coeff!!(a, 0.5)
    @test c !== a
    @test c == f1 + 0.5
    @test NC.additive_coeff(c) isa Float64
    @test a == f1 + 3
end

@testitem "NCAdd: constructor without filtering" setup = [Fermions] begin
    import NonCommutativeProducts as NC
    f1 = Fermion(1)
    k = NC.NCMul(1, [f1])
    x = NC.NCAdd(0, Dict(k => 0); filter_zeros=false, filter_scalars=false)
    @test length(x.dict) == 1
    @test only(values(x.dict)) == 0
    @test NC.NCAdd(0, Dict(k => 0)) == 0
    # the terms are always a Dict with the coefficient type as value type
    @test typeof(x.dict) == Dict{typeof(k),Int}
end

@testitem "NCAdd: zero terms are filtered on all paths" setup = [Fermions, Bosons] begin
    using VectorInterface
    import NonCommutativeProducts as NC
    import NonCommutativeProducts: add!!, add!, scale!!
    f1, f2 = Fermion(1), Fermion(2)
    b1 = Boson()

    function iscanonical(y, ref)
        y == ref && hash(y) == hash(ref) && !any(iszero, values(y.dict))
    end

    @test iscanonical(add!!(f1 + 1, -f1), 1)
    @test iscanonical(add!!(f1 + 1, f1, -1, 1), 1)
    @test iscanonical(add!!(f1 + 1, -f1 + 0), 1)
    @test iscanonical(add!!(f1 + 1, -1 * (f1 + 1)), 0)
    @test iscanonical(add!!(f1 + 1, f1, 1, 0), f1)
    @test iscanonical(add!!(f1 + f2, -1.0 * f1), f2) # widening coefficients
    @test iscanonical(add!!(1 + f2, NC.NCMul(0, [f1])), 1 + f2) # no new zero term
    @test iscanonical(add!!(f1 + 1, NC.NCMul(2, typeof(f1)[])), f1 + 3) # scalar products go to the coefficient
    @test iscanonical(add!(f1 + 1, 1, 1, 0), 1)
    @test iscanonical(NC.scale!(f1 + 1, 0), 0)
    @test iscanonical(scale!!(f1 + 1, 0), 0)
    @test iscanonical(scale!!(f1 + 1, 0.0), 0)
    @test iscanonical(VectorInterface.scale!!(f1 + 1, VectorInterface.Zero()), 0)
    @test iscanonical(VectorInterface.scale!(f1 + 1, f1 + f2, 0), 0)
    @test iscanonical(VectorInterface.add!(f1 + 1, f1 + 0, -1, 1), 1)
    @test iscanonical(VectorInterface.add!!(f1 + 1, f1 + 0, -1, 1), 1)
    @test iscanonical(f1 + 1 - f1, 1)
    @test iscanonical((f1 + 1) + (-f1 + 0), 1)
    @test iscanonical(0 * (f1 + 1), 0)
    # mixed species: the keys of the terms have different types
    @test iscanonical(add!!(1 * f1 + b1, -b1), f1 + 0)
    @test iscanonical(add!!(f1 + 1, b1, 1, 0), b1 + 0)
    @test iscanonical(add!!(1 * f1 + b1 + 1, -f1 - b1), 1)
    @test iscanonical(1 * f1 + b1 - b1, f1 + 0)
end

@testitem "NCAdd: unfilterable zero terms are kept by add!!" begin
    import NonCommutativeProducts as NC
    struct UnfilterableAtom end
    NC.isfilterable(::UnfilterableAtom) = false
    g = NC.NCMul(1, [UnfilterableAtom()])
    y = NC.add!!(g + 1, NC.NCMul(-1, [UnfilterableAtom()]))
    @test length(y.dict) == 1
    @test only(values(y.dict)) == 0
    y = NC.scale!!(g + 1, 0)
    @test length(y.dict) == 1
    @test iszero(y)
end

@testitem "NCAdd: UniformScaling + NCAdd copies" setup = [Fermions] begin
    using LinearAlgebra
    import NonCommutativeProducts: add!!
    f1, f2 = Fermion(1), Fermion(2)
    for op in (x -> I + x, x -> x + I, x -> 2I + x, x -> x - I)
        x = f1 + 1
        y = op(x)
        @test y !== x
        @test y.dict !== x.dict
        add!!(y, f2)
        @test x == f1 + 1
    end
end

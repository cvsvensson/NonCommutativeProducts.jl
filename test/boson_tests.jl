@testitem "Boson normal ordering" setup = [Bosons] begin
    import NonCommutativeProducts: NCMul, NCAdd, additive_coeff, prefactor
    using LinearAlgebra
    NonCommutativeProducts.enable_autosort!()
    b = Boson()

    @test b * b' == b' * b + 1
    @test b^2 * b' == b' * b^2 + 2 * b
    @test b * b'^2 == b'^2 * b + 2 * b'
    @test b^2 * b'^2 == b'^2 * b^2 + 4 * b' * b + 2

    # Compare with truncated matrices. Products of total degree d are exact on the first N - d Fock states.
    N = 16
    B = diagm(1 => sqrt.(1:N-1))
    mat(x::Number) = x * I(N)
    mat(x::Boson) = x.exp < 0 ? B^(-x.exp) : B'^x.exp
    mat(x::NCMul) = prefactor(x) * prod(mat, x.factors; init=Matrix(1.0I, N, N))
    mat(x::NCAdd) = additive_coeff(x) * I(N) + sum(v * mat(k) for (k, v) in pairs(x.dict); init=zeros(N, N))
    block = 1:N-8
    for n in 1:4, m in 1:4
        x = b^n * b'^m
        @test mat(x)[block, block] ≈ (B^n * B'^m)[block, block]
    end
    for word in ([b, b', b, b'], [b, b, b', b, b'], [b', b, b, b', b', b])
        x = prod(word)
        @test mat(x)[block, block] ≈ prod(mat, word)[block, block]
    end
end

@testitem "Boson Fock states" setup = [Bosons] begin
    using LinearAlgebra
    using VectorInterface: inner
    NonCommutativeProducts.enable_autosort!()
    b = Boson()
    ket(n) = State(n)
    bra(n) = State(n)'

    @test b' * ket(0) == ket(1)
    @test b * ket(0) == 0
    @test b * ket(1) == ket(0)
    @test b' * ket(1) == sqrt(2) * ket(2)
    @test b'^2 * ket(1) == sqrt(2) * sqrt(3) * ket(3)
    @test b^2 * ket(1) == 0
    @test bra(1) * b' == bra(0)
    @test bra(1) * b == sqrt(2) * bra(2)
    @test bra(0) * b' == 0
    @test (b' * ket(2))' == bra(2) * b

    @test inner(ket(2), ket(2)) == 1
    @test inner(ket(2), ket(3)) == 0
    for n in 0:5
        @test inner(ket(n), b' * b * ket(n)) ≈ n
        @test inner(b' * ket(n), b' * ket(n)) ≈ n + 1
        @test inner(ket(n), (b * b' - b' * b) * ket(n)) ≈ 1
    end

    # Matrix elements agree with truncated matrices, away from the truncation
    N = 16
    B = diagm(1 => sqrt.(1:N-1))
    for (x, X) in ((b^2 * b'^3 * b, B^2 * B'^3 * B), (b' * b^2 + b'^2, B' * B^2 + B'^2), (b * b' * b * b', B * B' * B * B'))
        @test [inner(ket(i), x * ket(j)) for i in 0:6, j in 0:6] ≈ X[1:7, 1:7]
    end

    # |n⟩⟨m| is an operator
    @test ket(1) * bra(0) * ket(0) == ket(1)
    @test ket(1) * bra(0) * ket(2) == 0
    @test bra(1) * (ket(1) * bra(0)) == bra(0)

    @test_throws ArgumentError b * bra(1)
    @test_throws ArgumentError ket(1) * b
    @test_throws ArgumentError ket(1) * ket(2)
    @test_throws ArgumentError bra(1) * bra(2)
end

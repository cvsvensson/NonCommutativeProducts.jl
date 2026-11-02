@testmodule GroupSpecies begin
    export GF, unhook
    import NonCommutativeProducts as NCP
    import NonCommutativeProducts: Swap, AddTerms
    # Two species of fermions, :a and :b, which anticommute with each other, and spins :s, which commute with both
    # and square to one. Species are ordered :a, :b, :s. H is :plain for no group hooks, :hooked for group_commutation,
    # and :ordered for group_commutation and group_isless.
    struct GF{S,H}
        label::Int
        creation::Bool
    end
    GF{S,H}(label) where {S,H} = GF{S,H}(label, false)
    Base.adjoint(x::GF{S,H}) where {S,H} = GF{S,H}(x.label, !x.creation)
    Base.show(io::IO, x::GF{S,H}) where {S,H} = print(io, S, x.creation ? "†" : "", x.label)
    NCP.@nc GF
    isfermion(S) = S !== :s
    rank(S) = findfirst(==(S), (:a, :b, :s))
    function NCP.mul_effect(x::GF{S,H}, y::GF{S,H}) where {S,H}
        isfermion(S) || return x.label > y.label ? Swap(1) : x.label == y.label ? 1 : nothing
        x == y && return 0
        (!x.creation, x.label) > (!y.creation, y.label) || return nothing
        x.label == y.label && return AddTerms((Swap(-1), 1))
        return Swap(-1)
    end
    function NCP.mul_effect(::GF{S1,H}, ::GF{S2,H}) where {S1,S2,H}
        rank(S1) < rank(S2) && return nothing
        return Swap(isfermion(S1) && isfermion(S2) ? -1 : 1)
    end
    NCP.group_commutation(::Type{GF{S1,H}}, ::Type{GF{S2,H}}) where {S1,S2,H} =
        H === :plain ? nothing : isfermion(S1) && isfermion(S2) ? :graded : :commute
    NCP.group_isless(::Type{GF{S1,:ordered}}, ::Type{GF{S2,:ordered}}) where {S1,S2} = rank(S1) < rank(S2)
    NCP.atom_parity(::GF{S}) where {S} = isfermion(S)
    unhook(x::GF{S}) where {S} = GF{S,:plain}(x.label, x.creation)
end

@testitem "Group hooks give the same sorted results" setup = [GroupSpecies] begin
    import NonCommutativeProducts as NCP
    import NonCommutativeProducts: NCMul, with_autosort, ncmap
    # The result of `f(H)` for species with the hooks H, mapped to the species without hooks
    plain(x) = with_autosort(() -> ncmap(unhook, x), false)
    function species(H)
        a, b, s = GF{:a,H}, GF{:b,H}, GF{:s,H}
        # same labels give extra terms (a a† = 1 - a† a, s s = 1), and a repeated fermion gives zero
        long = [s(2), b(1)', a(3), s(1), b(1), a(1)', a(3)', s(2), b(2), a(2), b(4)', a(5), s(3), b(3)', a(4)', s(4)]
        zero_product = [a(1), b(1), s(1), b(2), s(2), a(2), b(3), a(3), s(3), b(4), a(4), s(4), b(5), a(5), s(5), a(1)]
        h = sum(a(n)' * b(n) * s(n) * a(n + 1) + s(n + 1) * b(n)' * a(n) * s(n) for n in 1:2)
        return long, zero_product, h
    end
    long, zero_product, h = species(:plain)
    expected = sort(NCMul(1, Any[long...]))
    @test NCP.nterms(expected) > 1
    for H in (:hooked, :ordered)
        hlong, hzero, hh = species(H)
        for F in (Any, GF, Union{map(typeof, hlong)...})
            x = NCMul(2, F[hlong...])
            @test plain(sort(x)) == 2 * expected
            @test sort(sort(x)) == sort(x)
        end
        @test plain(sort(NCMul(1, Tuple(hlong)))) == expected
        @test iszero(sort(NCMul(1, Any[hzero...])))
        @test iszero(sort(NCMul(0, Any[hlong...])))
        # products of 16 factors
        @test plain(with_autosort(() -> hh^4, true)) == with_autosort(() -> h^4, true)
        @test plain(sort(with_autosort(() -> hh^4, false))) == sort(with_autosort(() -> h^4, false))
    end
end

@testitem "Graded group commutation" setup = [GroupSpecies] begin
    import NonCommutativeProducts: NCMul
    a, b, s = GF{:a,:hooked}, GF{:b,:hooked}, GF{:s,:hooked}
    as = [a(n) for n in 1:15]
    # b1 moves past 15 anticommuting fermions, and s1 past 15 commuting ones
    @test sort(NCMul(1, Any[b(1), as...])) == NCMul(-1, Any[as..., b(1)])
    @test sort(NCMul(1, Any[s(1), as...])) == NCMul(1, Any[as..., s(1)])
    # b1 b2 moves past 14 fermions as a pair, and s1 moves past both kinds
    @test sort(NCMul(1, Any[b(1), b(2), s(1), as[1:14]...])) == NCMul(1, Any[as[1:14]..., b(1), b(2), s(1)])
    @test sort(NCMul(1, Any[b(1), as..., a(16)])) == NCMul(1, Any[as..., a(16), b(1)])
end

@testitem "Groups without declared commutation keep their order" begin
    import NonCommutativeProducts as NCP
    import NonCommutativeProducts: NCMul, Swap
    # Commuting atoms of two types that sort by label across types, so the types interleave
    struct LeftLabel
        n::Int
    end
    struct RightLabel
        n::Int
    end
    NCP.@nc LeftLabel RightLabel
    const Label = Union{LeftLabel,RightLabel}
    NCP.mul_effect(x::Label, y::Label) = x.n > y.n ? Swap(1) : nothing
    labels = [4, 3, 2, 1, 8, 7, 6, 5, 12, 11, 10, 9, 16, 15, 14, 13]
    atom(n) = isodd(n) ? LeftLabel(n) : RightLabel(n)
    @test sort(NCMul(1, Any[atom(n) for n in labels])) == NCMul(1, Any[atom(n) for n in 1:16])
end

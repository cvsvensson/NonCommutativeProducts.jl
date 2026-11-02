@testmodule GroupSpecies begin
    export GF, unhook
    import NonCommutativeProducts as NCP
    import NonCommutativeProducts: Swap, AddTerms
    # Two species of fermions, :a and :b, which anticommute with each other, and spins :s, which commute with both
    # and square to one. Species are ordered :a, :b, :s. H is :plain for no group hooks and :hooked for
    # group_commutation.
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
    NCP.atom_parity(::GF{S}) where {S} = isfermion(S)
    unhook(x::GF{S}) where {S} = GF{S,:plain}(x.label, x.creation)
end

@testitem "Group hooks give the same sorted results" setup = [GroupSpecies] begin
    import NonCommutativeProducts as NCP
    import NonCommutativeProducts: NCMul, with_autosort, ncmap
    plain(x) = with_autosort(() -> ncmap(unhook, x), false)
    # long enough to sort by group; same labels give extra terms (a a† = 1 - a† a, s s = 1)
    n = NCP._group_sort_min_length
    species(H) = (GF{:a,H}, GF{:b,H}, GF{:s,H})
    function long(H)
        a, b, s = species(H)
        block(o) = [s(o + 2), b(o + 1)', a(o + 3), s(o + 1), b(o + 1), a(o + 1)', a(o + 3)', s(o + 2), b(o + 2),
            a(o + 2), b(o + 4)', a(o + 5), s(o + 3), b(o + 3)', a(o + 4)', s(o + 4)]
        return reduce(vcat, block(10j) for j in 0:cld(n, 16)-1)
    end
    expected = sort(NCMul(1, Any[long(:plain)...]))
    @test NCP.nterms(expected) > 1
    for F in (Any, GF)
        x = NCMul(2, F[long(:hooked)...])
        @test plain(sort(x)) == 2 * expected
        @test sort(sort(x)) == sort(x)
    end
    a, b, s = species(:hooked)
    @test iszero(sort(NCMul(1, Any[long(:hooked)..., a(1), a(1)])))
    # products of 16 factors from powers of a sum
    h((A, B, S)) = sum(A(k)' * B(k) * S(k) * A(k + 1) + S(k + 1) * B(k)' * A(k) * S(k) for k in 1:2)
    @test plain(with_autosort(() -> h(species(:hooked))^4, true)) == with_autosort(() -> h(species(:plain))^4, true)
end

@testitem "Graded group commutation" setup = [GroupSpecies] begin
    import NonCommutativeProducts as NCP
    import NonCommutativeProducts: NCMul
    a, b, s = GF{:a,:hooked}, GF{:b,:hooked}, GF{:s,:hooked}
    m = NCP._group_sort_min_length
    as = [a(n) for n in 1:m]
    # b1 moves past an even number of anticommuting fermions, and past an odd number, and s1 past commuting ones
    @test sort(NCMul(1, Any[b(1), as...])) == NCMul(iseven(m) ? 1 : -1, Any[as..., b(1)])
    @test sort(NCMul(1, Any[b(1), as..., a(m + 1)])) == NCMul(iseven(m) ? -1 : 1, Any[as..., a(m + 1), b(1)])
    @test sort(NCMul(1, Any[s(1), as...])) == NCMul(1, Any[as..., s(1)])
    # b1 b2 moves past the fermions as a pair, and s1 moves past both kinds
    @test sort(NCMul(1, Any[b(1), b(2), s(1), as...])) == NCMul(1, Any[as..., b(1), b(2), s(1)])
end

@testitem "Default hooks keep result types" setup = [GroupSpecies] begin
    import NonCommutativeProducts as NCP
    import NonCommutativeProducts: NCMul
    x = NCMul(1, Any[GF{:a,:plain}(n) for n in 2 * NCP._group_sort_min_length:-1:1])
    @test typeof(sort(x)) == typeof(sort(NCMul(1, Any[GF{:a,:plain}(1), GF{:a,:plain}(2)])))
end

@testitem "Inconsistent group hooks throw" begin
    import NonCommutativeProducts as NCP
    import NonCommutativeProducts: NCMul, Swap
    struct HookX
        n::Int
    end
    struct HookY
        n::Int
    end
    NCP.@nc HookX HookY
    NCP.mul_effect(x::T, y::T) where {T<:Union{HookX,HookY}} = x.n > y.n ? Swap(1) : nothing
    NCP.mul_effect(::HookX, ::HookY) = nothing
    NCP.mul_effect(::HookY, ::HookX) = Swap(1)
    NCP.atom_parity(::Union{HookX,HookY}) = true
    m = 2 * NCP._group_sort_min_length
    atom(n) = isodd(n) ? HookY(n) : HookX(n)
    x = NCMul(1, Any[atom(n) for n in 1:m])
    # the commutation declared for (HookX, HookY) and for (HookY, HookX)
    const declared = Ref{Tuple{Any,Any}}((nothing, :graded))
    NCP.group_commutation(::Type{HookX}, ::Type{HookY}) = declared[][1]
    NCP.group_commutation(::Type{HookY}, ::Type{HookX}) = declared[][2]
    @test_throws "not symmetric" sort(x)
    # declared graded, but mul_effect swaps odd atoms with Swap(1)
    declared[] = (:graded, :graded)
    @test_throws "must be nothing or Swap(-1)" sort(x)
    declared[] = (:commute, :commute)
    @test sort(x) == NCMul(1, Any[atom(n) for n in [2:2:m; 1:2:m]])
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

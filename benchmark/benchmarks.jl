using BenchmarkTools
using NonCommutativeProducts
using NonCommutativeProducts: with_autosort, NCMul, NCAdd
using LinearAlgebra: norm
using Random
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1
const SUITE = BenchmarkGroup()
Random.seed!(1)
const VI = NonCommutativeProducts.VectorInterface

include("../test/rules/fermions.jl")
include("../test/rules/bosons.jl")
NonCommutativeProducts.@commutative Boson Fermion
NonCommutativeProducts.enable_autosort!()

# Building sums: each `+` of two sums copies the terms of the left one
SUITE["construct"] = BenchmarkGroup()
SUITE["construct"]["sum"] = @benchmarkable sum(Fermion(n)' * Fermion(n) + Fermion(n) * Fermion(n + 1)' for n in 1:1000)

# Products of sums and atoms, sorted as they are formed
SUITE["multiply"] = BenchmarkGroup()
const h = sum(Fermion(n)' * Fermion(n) + Fermion(n) * Fermion(n + 1)' for n in 1:10)
SUITE["multiply"]["sum_cube"] = @benchmarkable $h^3

# The same sum, for fermions whose `mul_effect` returns `Swap` and `AddTerms` instead of forming new expressions
struct SwapFermion
    label::Int
    creation::Bool
end
SwapFermion(k) = SwapFermion(k, false)
Base.adjoint(x::SwapFermion) = SwapFermion(x.label, !x.creation)
Base.hash(x::SwapFermion, h::UInt) = hash(x.label, hash(x.creation, hash(0x452821e638d01377, h)))
NonCommutativeProducts.@nc SwapFermion
function NonCommutativeProducts.mul_effect(a::SwapFermion, b::SwapFermion)
    a == b && return 0
    (!a.creation, a.label) > (!b.creation, b.label) || return nothing
    a.label == b.label && return NonCommutativeProducts.AddTerms((NonCommutativeProducts.Swap(-1), 1))
    return NonCommutativeProducts.Swap(-1)
end
const h_swap = sum(SwapFermion(n)' * SwapFermion(n) + SwapFermion(n) * SwapFermion(n + 1)' for n in 1:10)
SUITE["multiply"]["sum_cube_swap"] = @benchmarkable $h_swap^3

labels = shuffle(1:10)
SUITE["multiply"]["deep_product"] = @benchmarkable prod(Fermion(l) for l in $labels) * prod(Fermion(l)' for l in $labels)

# Normal ordering of a single bosonic mode, where every b b† gives an AddTerms with a Swap
SUITE["multiply"]["boson_power"] = @benchmarkable (Boson() + Boson()')^12

# Products of bosons and fermions have factors in a Vector{Any}. Boson is a single mode, and Boson() * Boson()' is the
# only boson pair that adds a term when reordered (b b† = b† b + 1). The Boson in the middle of the last term is moved
# to the left of the Fermion, which is a plain swap.
mixed_term(n) = Fermion(n)' * Fermion(n) + Fermion(n) * Fermion(n + 1)' + Boson()' * Fermion(n) * Boson()
const h_mixed = sum(mixed_term(n) for n in 1:8)
SUITE["multiply"]["mixed_sum_cube"] = @benchmarkable $h_mixed^3

# A single product with lots of reordering but few additions: fermions with distinct labels only swap, and bosons
# only swap past fermions, so the extra terms come from the three b b† pairs.
mixed_ops = shuffle([Fermion.(1:6); adjoint.(Fermion.(7:12)); fill(Boson(), 3); fill(Boson()', 3)])
SUITE["multiply"]["mixed_deep_product"] = @benchmarkable prod($mixed_ops)

# Without autosort, the default: products are added up unsorted, and sorted only by an explicit `sort`
SUITE["unsorted"] = BenchmarkGroup()
SUITE["unsorted"]["sum_cube"] = @benchmarkable with_autosort(() -> $h^3, false)
const h3_unsorted = with_autosort(() -> h^3, false)
SUITE["unsorted"]["sort"] = @benchmarkable sort($h3_unsorted)

# The commutator of an interacting chain with the position operator. The density and interaction terms cancel, and
# the current terms remain.
SUITE["operators"] = BenchmarkGroup()
const H_int = sum(Fermion(n)' * Fermion(n + 1) + Fermion(n + 1)' * Fermion(n) +
    0.5 * Fermion(n)' * Fermion(n) * Fermion(n + 1)' * Fermion(n + 1) for n in 1:30)
const X = sum(n * Fermion(n)' * Fermion(n) for n in 1:30)
SUITE["operators"]["commutator"] = @benchmarkable $H_int * $X - $X * $H_int

# Expressions as vectors: a tight-binding chain, with states as sums of kets and the Hamiltonian as a sum of |n⟩⟨m|.
# Applying it to a state reduces each ⟨m|k⟩ to a scalar.
struct Ket
    n::Int
    bra::Bool
end
Ket(n) = Ket(n, false)
Base.adjoint(k::Ket) = Ket(k.n, !k.bra)
Base.hash(k::Ket, h::UInt) = hash(k.n, hash(k.bra, hash(0x13198a2e03707344, h)))
NonCommutativeProducts.@nc Ket
NonCommutativeProducts.mul_effect(a::Ket, b::Ket) = a.bra && !b.bra ? Int(a.n == b.n) : nothing

# A bicyclic grade, so that zero sums of kets are zero kets
struct Bicyclic
    m::Int
    n::Int
end
NonCommutativeProducts.grade(k::Ket) = k.bra ? Bicyclic(0, 1) : Bicyclic(1, 0)
function NonCommutativeProducts.mul_grade(g::Bicyclic, h::Bicyclic)
    c = max(g.n, h.m)
    r = Bicyclic(g.m - g.n + c, h.n - h.m + c)
    return r == Bicyclic(0, 0) ? nothing : r
end
NonCommutativeProducts.adjoint_grade(g::Bicyclic) = Bicyclic(g.n, g.m)

const N_chain = 50
const H_chain = sum(Ket(n) * Ket(n + 1)' + Ket(n + 1) * Ket(n)' for n in 1:(N_chain-1)) + sum(0.1n * Ket(n) * Ket(n)' for n in 1:N_chain)
const x_chain = sum(randn() * Ket(n) for n in 1:N_chain)
const y_chain = sum(randn() * Ket(n) for n in 1:N_chain)

# An orthonormal Krylov basis, built with the VectorInterface operations that iterative solvers use
function arnoldi(H, x, n)
    basis = [VI.scale(x, inv(norm(x)))]
    for _ in 2:n
        w = H * basis[end]
        for v in basis
            w = VI.add!!(w, v, -VI.inner(v, w), VI.One())
        end
        push!(basis, VI.scale!!(w, inv(norm(w))))
    end
    return basis
end

SUITE["vector"] = BenchmarkGroup()
SUITE["vector"]["apply"] = @benchmarkable $H_chain * $x_chain
SUITE["vector"]["inner"] = @benchmarkable VI.inner($x_chain, $y_chain)
SUITE["vector"]["arnoldi"] = @benchmarkable arnoldi($H_chain, $x_chain, 10)

# The same chain with many more sites, where the work per pair of terms dominates
const N_large = 600
const H_large = sum(Ket(n) * Ket(n + 1)' + Ket(n + 1) * Ket(n)' for n in 1:(N_large-1)) + sum(0.1n * Ket(n) * Ket(n)' for n in 1:N_large)
const x_large = sum(randn() * Ket(n) for n in 1:N_large)
SUITE["vector"]["apply_large"] = @benchmarkable $H_large * $x_large

# Kets with a contraction key, so that `inner` pairs only terms with the same label
struct KeyedKet
    n::Int
    bra::Bool
end
KeyedKet(n) = KeyedKet(n, false)
Base.adjoint(k::KeyedKet) = KeyedKet(k.n, !k.bra)
Base.hash(k::KeyedKet, h::UInt) = hash(k.n, hash(k.bra, hash(0xa4093822299f31d0, h)))
NonCommutativeProducts.@nc KeyedKet
NonCommutativeProducts.mul_effect(a::KeyedKet, b::KeyedKet) = a.bra && !b.bra ? Int(a.n == b.n) : nothing
NonCommutativeProducts.contraction_key(k::KeyedKet) = k.n

const N_keyed = 4096
const x_keyed = sum(randn() * KeyedKet(n) for n in 1:N_keyed)
const y_keyed = sum(randn() * KeyedKet(n) for n in 1:N_keyed)
SUITE["vector"]["inner_keyed"] = @benchmarkable VI.inner($x_keyed, $y_keyed)

# |n⟩⟨m| applied to a sum of kets is the coefficient of |m⟩ times |n⟩, found by a scan over the terms of the sum.
# Other products go through the default method.
_isket(fs) = length(fs) == 1 && !only(fs).bra
function NonCommutativeProducts.mul_terms!!(c, x::NCMul{<:Any,KeyedKet}, b::NCAdd)
    fs = x.factors
    if length(fs) != 2 || fs[1].bra || !fs[2].bra || !all(_isket ∘ first, NonCommutativeProducts.termpairs(b))
        return invoke(NonCommutativeProducts.mul_terms!!, Tuple{Any,NCMul,NCAdd}, c, x, b)
    end
    n, m = fs[1].n, fs[2].n
    for (kets, coeff) in NonCommutativeProducts.termpairs(b)
        if only(kets).n == m
            c = NonCommutativeProducts.add!!(c, NonCommutativeProducts.coefficient(x) * coeff * KeyedKet(n))
        end
    end
    return c
end

const H_keyed = sum(KeyedKet(n) * KeyedKet(n + 1)' + KeyedKet(n + 1) * KeyedKet(n)' for n in 1:(N_large-1)) +
    sum(0.1n * KeyedKet(n) * KeyedKet(n)' for n in 1:N_large)
const x_keyed_large = sum(randn() * KeyedKet(n) for n in 1:N_large)
SUITE["vector"]["mul_hook"] = @benchmarkable $H_keyed * $x_keyed_large

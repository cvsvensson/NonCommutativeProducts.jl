using NonCommutativeProducts
import NonCommutativeProducts: @nc, AddTerms, Swap, NCMul
struct Boson
    exp::Int
end
Base.adjoint(x::Boson) = Boson(-x.exp)
Boson() = Boson(-1)
# A session-independent hash, see the comment on the hash of Fermion
Base.hash(x::Boson, h::UInt) = hash(x.exp, hash(0x2e8d4b7a1c9f3e05, h))
Base.show(io::IO, x::Boson) = print(io, "b", x.exp > 0 ? "†" : "", abs(x.exp) > 1 ? "^($(x.exp))" : "")
@nc Boson

function _normal_term(n, m, k)
    # checked, since the product overflows Int already for n = m - 1 = 19
    coeff = Base.checked_mul(factorial(k), binomial(n, k), binomial(m, k))
    factors = [Boson(e) for e in (m - k, k - n) if e != 0]
    return isempty(factors) ? coeff : NCMul(coeff, factors)
end
function NonCommutativeProducts.mul_effect(a::Boson, b::Boson)
    sign(a.exp) == sign(b.exp) && return Boson(a.exp + b.exp)
    if a.exp < 0 && b.exp > 0
        # Normal ordering: b^n b'^m = Σ_k k! binomial(n,k) binomial(m,k) b'^(m-k) b^(n-k), with the k = 0 term a swap
        n, m = -a.exp, b.exp
        return AddTerms((Swap(1), (_normal_term(n, m, k) for k in 1:min(n, m))...))
    else
        return nothing
    end
end

# A Fock state |occ⟩ of the single mode, or the bra ⟨occ| when adj is true
struct BState
    occ::Int
    adj::Bool
end
State(occ) = BState(occ, false)
Base.adjoint(s::BState) = BState(s.occ, !s.adj)
Base.hash(s::BState, h::UInt) = hash(s.occ, hash(s.adj, hash(0x7c1e5a9d3b2f8064, h)))
Base.show(io::IO, s::BState) = print(io, s.adj ? "⟨$(s.occ)|" : "|$(s.occ)⟩")

# b'^m|n⟩ = √((n+m)!/n!) |n+m⟩ and b^m|n⟩ = √(n!/(n-m)!) |n-m⟩. On a bra, b' lowers and b raises the occupation.
function _apply(s::BState, change)
    newocc = s.occ + change
    newocc < 0 && return 0
    coeff = prod(sqrt, min(s.occ, newocc)+1:max(s.occ, newocc); init=1.0)
    return NCMul(coeff, [BState(newocc, s.adj)])
end
function NonCommutativeProducts.mul_effect(b::Boson, s::BState)
    s.adj && throw(ArgumentError("Cannot multiply Boson by a bra from the left"))
    return _apply(s, b.exp)
end
function NonCommutativeProducts.mul_effect(s::BState, b::Boson)
    s.adj || throw(ArgumentError("Cannot multiply a ket by Boson from the right"))
    return _apply(s, -b.exp)
end
function NonCommutativeProducts.mul_effect(s1::BState, s2::BState)
    s1.adj && !s2.adj && return Int(s1.occ == s2.occ) # ⟨n|m⟩
    !s1.adj && s2.adj && return nothing # |n⟩⟨m| is an operator
    throw(ArgumentError("Cannot multiply two kets or two bras of the same mode"))
end
NonCommutativeProducts.@nc BState Boson

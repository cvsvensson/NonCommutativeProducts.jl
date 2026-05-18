using NonCommutativeProducts
import NonCommutativeProducts: @nc

"""
    Pauli(label, idx)

Pauli matrix operator for spin-1/2 systems.
- `label`: identifies the site/particle
- `idx`: 1 for σx, 2 for σy, 3 for σz
"""
struct Pauli{L}
    label::L
    idx::Int
    Pauli(label::L, idx::Int=1) where L = new{L}(label, idx)
end

Base.adjoint(x::Pauli) = x  # Hermitian
Pauli(idx::Int) = Pauli(0, idx)

function Base.show(io::IO, x::Pauli)
    syms = ("x", "y", "z")
    base = "σ" * syms[x.idx]
    print(io, x.label == 0 ? base : "$base[$(x.label)]")
end

@nc Pauli

function should_swap(a::Pauli, b::Pauli)
    a.label != b.label ? a.label > b.label : a.idx > b.idx
end

function NonCommutativeProducts.mul_effect(a::Pauli, b::Pauli)
    if a.label != b.label
        # Different sites commute
        return should_swap(a, b) ? b * a : nothing
    elseif a.idx == b.idx
        # σᵢ² = 1
        return 1
    elseif should_swap(a, b)
        # σᵢ σⱼ = -σⱼ σᵢ for i ≠ j (anti-commute)
        return -b * a 
    else
        return nothing
    end
end
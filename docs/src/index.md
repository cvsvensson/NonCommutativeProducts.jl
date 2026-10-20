# NonCommutativeProducts.jl

**NonCommutativeProducts.jl** is a Julia package for sorting non-commuting objects, such as operators in quantum mechanics. Users must specify custom commutation relations and sorting orders, as there are no inbuilt ones in this package.

## Installation

```julia
using Pkg
Pkg.add("NonCommutativeProducts")
```

## How to do it

With a list of types `Ts` that represents your non-commuting objects, call [`@nc Ts...`](@ref NonCommutativeProducts.@nc) which defines addition and multiplication for these types. In order to sort them, define [`mul_effect(a::Tj, b::Tk)`](@ref NonCommutativeProducts.mul_effect) for each pair of types to define the behaviour of `a*b`.

The function `mul_effect(a::Tj, b::Tk)` can have return values
* `nothing`: Keeps `a*b`. This return value is important as the sorting only terminates when this is the return value for each neighbouring pair product. Without it, sorting never terminates.
* `λ::Number`: replaces `a*b` by `λ`.
* `x::T`: replaces `a*b` by `x`.
* Sums and products of such terms, such as `b*a + 1`.

No names are exported from the package, so import the names you need.

The example below shows how to define and sort fermions. The tests of this package contain a similar example with Majorana fermions, and another test shows how one can collect powers.

## Example: Fermions

Fermions satisfy
```math
\begin{align*}
\{c_i^\dagger,c_j\}&=\delta_{ij}\\
\{c_i,c_j\}&=0\\
\{c_i^\dagger,c_j^\dagger\}&=0.
\end{align*}
```

We'll sort them in normal order, i.e. all creation operators appear before annihilation operators. First, define a struct representing a fermion.
```jldoctest fermions; output = false
struct Fermion
    label::Int
    dagger::Bool
end
Base.adjoint(x::Fermion) = Fermion(x.label, !x.dagger)
Fermion(k) = Fermion(k, false)
Base.show(io::IO, x::Fermion) = print(io, "c", x.dagger ? "†" : "", "[", x.label, "]")

# output

```
Then, hook up the struct to the package to let it handle the arithmetic. Load the package and import the functions we are going to use.
```jldoctest fermions; output = false
using NonCommutativeProducts
import NonCommutativeProducts: @nc, mul_effect
@nc Fermion

# output

```
Now one can add and multiply these fermions,
```jldoctest fermions; setup = :(NonCommutativeProducts.disable_autosort!())
julia> Fermion(1)'*Fermion(1) + 1
1I + c†[1]*c[1]
```
but they can't be sorted. To sort them in normal order and by label, define
```jldoctest fermions; output = false, filter = r"\d+ methods?"
function mul_effect(a::Fermion, b::Fermion)
    # If the fermion is multiplied with itself, replace it by zero.
    (a.dagger, a.label) == (b.dagger, b.label) && return 0
    # If already sorted, return nothing
    (!a.dagger, a.label) < (!b.dagger, b.label) && return nothing
    # Non-trivial anti-commutation relation
    (!a.dagger, a.label) == (b.dagger, b.label) && return -b*a + 1
    # Trivial anti-commutation relation
    return -b*a
end

# output

mul_effect (generic function with 1 method)
```
Now we can sort expressions involving fermions while respecting the commutation relations.
```jldoctest fermions
julia> Fermion(1)'*Fermion(1) |> sort
c†[1]*c[1]

julia> Fermion(1)*Fermion(1)' |> sort
1I - c†[1]*c[1]
```
In order to automatically sort them on each multiplication, call [`enable_autosort!`](@ref NonCommutativeProducts.enable_autosort!):
```jldoctest fermions
julia> NonCommutativeProducts.enable_autosort!();

julia> prod(Fermion(n) + Fermion(n)' for n in 1:4)
Sum with 16 terms:
-c†[1]*c†[2]*c†[4]*c[3]-c†[1]*c†[3]*c[2]*c[4]-c†[2]*c[1]*c[3]*c[4] + ...
```
`enable_autosort!` sets the global default. You can override it locally in a scope with [`NonCommutativeProducts.with_autosort(false) do ... end`](@ref NonCommutativeProducts.with_autosort). This temporary override does not change the global default. When the function `mul_effect` is called from within this package, autosort is always locally disabled to avoid infinite recursion. Use [`NonCommutativeProducts.disable_autosort!()`](@ref NonCommutativeProducts.disable_autosort!) to turn autosort off globally again.
```jldoctest fermions
julia> NonCommutativeProducts.disable_autosort!();
```

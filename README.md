# NonCommutativeProducts.jl

<!--  [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cvsvensson.github.io/NonCommutativeProducts.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cvsvensson.github.io/NonCommutativeProducts.jl/dev/) -->
[![Build Status](https://github.com/cvsvensson/NonCommutativeProducts.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cvsvensson/NonCommutativeProducts.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/cvsvensson/NonCommutativeProducts.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cvsvensson/NonCommutativeProducts.jl)

**NonCommutativeProducts.jl** is a Julia package for sorting non-commuting objects, such as operators in quantum mechanics. Users must specify custom commutation relations and sorting orders, as there are no inbuilt ones in this package. 

## How to do it
With a list of types `Ts` that represents your non-commuting objects, call `@nc Ts...` which defines addition and multiplication for these types. In order to sort them, define `mul_effect(a::Tj, b::Tk)` for each pair of types to define the behaviour of `a*b`.

The function `mul_effect(a::Tj, b::Tk)` can have return values 
* `nothing`: Keeps `a*b`. This return value is important as the sorting only terminates when this is the return value for each neighbouring pair product. If you don't have this, you'll get stuck in an infinite loop.
* `λ::Number`: replaces `a*b` by `λ`.
* `x::T`: replaces `a*b` by `x`
* Sums and products of such terms, such as `b*a + 1`.

No names are currently exported from the package, so you'll need to import the names you need.

Let's look at an example that shows how to define and sort fermions. The tests of this package contains a similar example with majorana fermions, and another test shows how one can collect powers.

## Example: Fermions

Let's see how to define fermions which satisfy
```math
\begin{align*}
\{c_i^\dagger,c_j\}&=\delta_{ij}\\ 
\{c_i,c_j\}&=0\\ 
\{c_i^\dagger,c_j^\dagger\}&=0.
\end{align*}
```

We'll sort them in normal order i.e. all creation operators appear before annihilation operators. First let's define a struct representing a fermion.
```julia
struct Fermion
	label::Int
    dagger::Bool
end
Base.adjoint(x::Fermion) = Fermion(x.label, !x.dagger)
Fermion(k) = Fermion(k, false)
Base.show(io::IO, x::Fermion) = print(io, "c", x.dagger ? "†" : "", "[", x.label, "]")
```
Then, we need to hook up our struct to the package to let it handle the arithmetic. Let's load the package and import the functions we are gonna use.
```julia
using NonCommutativeProducts
import NonCommutativeProducts: @nc, mul_effect
@nc Fermion
```
Now one can add and multiply these fermions,
```julia
Fermion(1)'*Fermion(1) + 1
#1I + c†[1]*c[1]
```
but they can't be sorted. To sort them in normal order and by label, we define 
```julia
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
```
Now we can sort expression involving fermions while respecting the commutation relations.
```julia
Fermion(1)'*Fermion(1) |> sort
#c†[1]*c[1]
Fermion(1)*Fermion(1)' |> sort
#1I - c†[1]*c[1]
```
In order to automatically sort them on each multiplication, we can call `enable_autosort!`:
```julia
NonCommutativeProducts.enable_autosort!()
prod(Fermion(n) + Fermion(n)' for n in 1:4) 
#=Sum with 16 terms: 
 -c†[1]*c†[2]*c†[4]*c[3] + c†[1]*c[2]*c[3]*c[4] + c†[1]*c†[3]*c†[4]*c[2] + ...=#
```
`enable_autosort!` sets the global default. You can override it locally in a scope with `Base.ScopedValues.with(NonCommutativeProducts._autosort => false) do ... end`. This temporary override does not change the global default. When the function `mul_effect` is called from within this package, autosort is always locally disabled to avoid infinite recursion.

## Performance tips

This package is flexible, but not very efficient. Sorting is done via bubble sort, which is convenient for this use case since it is based on repeatedly swapping adjacent elements where commutation relations can be used. But it does not scale well with the length of the list, so it won't perform well for products of many elements.

```julia
# Example timing: 85.840 ms (1132411 allocations: 82.81 MiB)
op = prod(Fermion(n) + Fermion(n)' + 1 for n in 1:10)
#= Sum with 59049 terms: 
1I-c†[1]*c†[2]*c†[6]*c[5]*c[7]*c[10]+c†[8]*c[2]*c[3]*c[5]*c[6]*c[9]*c[10]+c†[1]*c†[4]*c†[6]*c†[8]*c†[10]*c[2]*c[5]*c[9] + ...=#
```


To cut down on allocations, one can use `NonCommutativeProducts.add!!` which tries to perform addition in place, but widens the type if not possible.
```julia
op = 0
for n in 1:100
    op += Fermion(n)'*Fermion(n)
end
# Example timing: 83.000 μs (2280 allocations: 403.52 KiB)

op2 = zero(op)
for n in 1:100
    op2 = NonCommutativeProducts.add!!(op2, Fermion(n)'*Fermion(n))
end
# Example timing: 36.200 μs (1708 allocations: 109.61 KiB)
op == op2 #true
```

Mixing different symbolic types together is slower than using a single type.

Disable autosort to improve performance if you don't want automatic sorting.

Finally, instead of using ordinary addition and multiplication in mul_effect, 
you can use `Swap` and `AddTerms` which are more efficient. 
Their behaviour is:
* `Swap(λ::Number)`: Replaces `a*b` by `λ*b*a`. 
* `AddTerms(terms)`: `a*b` should be replaced by a sum of terms. `terms` should be an iterable such as a vector or a tuple, and the elements can be of the other allowed return types.

For example, the `mul_effect` for fermions can be implemented as
```julia
import NonCommutativeProducts: Swap, AddTerms
function mul_effect(a::Fermion, b::Fermion)
    (a.dagger, a.label) == (b.dagger, b.label) && return 0 
    (!a.dagger, a.label) < (!b.dagger, b.label) && return nothing
    (!a.dagger, a.label) == (b.dagger, b.label) && return AddTerms((Swap(-1), 1))
    return Swap(-1)
end
```

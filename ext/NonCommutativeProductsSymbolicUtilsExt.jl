module NonCommutativeProductsSymbolicUtilsExt

using NonCommutativeProducts
import SymbolicUtils

const NC = NonCommutativeProducts

# Substitutes in the coefficients only; atoms are left as they are.
SymbolicUtils.substitute(x::NC.MulAdd, dict; kwargs...) =
    NC.map_coefficients((_, c) -> SymbolicUtils.substitute(c, dict; kwargs...), x)

# Symbolic constants have the bound of their value, and other symbolic coefficients have no bound.
NC._unwrap_symbolic(c) = SymbolicUtils.unwrap(c)
NC.coefficient_bound(c::SymbolicUtils.BasicSymbolic) =
    SymbolicUtils.isconst(c) ? NC.coefficient_bound(SymbolicUtils.unwrap_const(c)) : nothing

end

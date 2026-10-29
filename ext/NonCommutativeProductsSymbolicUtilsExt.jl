module NonCommutativeProductsSymbolicUtilsExt

using NonCommutativeProducts
import SymbolicUtils

const NC = NonCommutativeProducts

# Substitutes in the coefficients only; atoms are left as they are.
SymbolicUtils.substitute(x::NC.MulAdd, dict; kwargs...) =
    NC.map_coefficients((_, c) -> SymbolicUtils.substitute(c, dict; kwargs...), x)

end

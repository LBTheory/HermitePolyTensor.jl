#------------------------------------------------------------------------------#
#                                   Package                                    #
#------------------------------------------------------------------------------#

"""
# 1. Description

    module HermitePolyTensor

Symbolic Hermite Polynomial Tensors — or N-dimensional Hermite polynomials [1] —
in one to three Euclidean space dimensions.

# 2. Usage

```julia-repl
julia> import HermitePolyTensor

```

# 3. References

[1]: H. Grad, “Note on N-dimensional  Hermite  polynomials,”  Communications  on
Pure and Applied Mathematics, vol. 2, no. 4, pp. 325–330, 1949.
"""
module HermitePolyTensor


#------------------------------------------------------------------------------#
#                                   Imports                                    #
#------------------------------------------------------------------------------#

using Reexport


#------------------------------------------------------------------------------#
#                                   Includes                                   #
#------------------------------------------------------------------------------#

include("version.jl")

include("HPT.jl")
@reexport using .HPT


end # module

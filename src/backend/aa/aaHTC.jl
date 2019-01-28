#------------------------------------------------------------------------------#
#                                    Module                                    #
#------------------------------------------------------------------------------#

"""
# Description

    module aaHTC

AbstractAlgebra backend (aa) module to the Hermite Tensor Polynomial package.

# Usage

This module is `include`d by the `HermitePolyTensor` package and is not  usually
used in isolation, or directly; rather through the  interfaces  defined  on  the
`HPT` submodule, which is loaded automatically by the package:

```julia-repl
julia> using HermitePolyTensor

julia> typeof(HPT.aaHTC)
Module

```
"""
module aaHTC


#------------------------------------------------------------------------------#
#                                   Imports                                    #
#------------------------------------------------------------------------------#

import AbstractAlgebra


#------------------------------------------------------------------------------#
#                                   Exports                                    #
#------------------------------------------------------------------------------#

# Auxiliary Functions
export SVHP 

# Tensor Component Functions
export HTC


#------------------------------------------------------------------------------#
#                             Auxiliary Functions                              #
#------------------------------------------------------------------------------#

"""
# Description

    SVHP(
        n::Int; ax::Int = 1, st::Bool = true
    )::AbstractAlgebra.Generic.MPoly{Rational{BigInt}}

Returns a symbolic `n`-th order single-variable Hermite polynomial in  the  axis
`ax`,  axis  1,  2,  and  3  corresponding  to  the  x,  y,  and  z  directions,
respectively. A "probabilist's" polynomial is returned if `st`  is  `true`  (the
default), otherwise, a "physicist's".

This function is not used directly; rather through the interfaces defined on the
`HPT` submodule, which is loaded automatically by the package, by specifying the
AbstractAlgebra ("aa") backend:

# Usage

```julia-repl
julia> using HermitePolyTensor

julia> HPT.SVHP(3, be = "aa")
x^3-3//1*x

julia> HPT.SVHP(4, ax = 2, be = "aa")
y^4-6//1*y^2+3//1

julia> HPT.aaHTC.SVHP(5, ax = 3, st = false) # Direct access
32//1*z^5-160//1*z^3+120//1*z

```
"""
function SVHP(
    n::Int; ax::Int = 1, st::Bool = true
)::AbstractAlgebra.Generic.MPoly{Rational{BigInt}}
    # Validations
    if n < 0
        throw(DomainError("n = $n outside the valid domain [0, ∞)"))
    elseif ax <= 0 || ax >= 4
        throw(DomainError("ax = $ax outside the valid domain [1, 3]"))
    end
    # Execution
    R, (x, y, z) = AbstractAlgebra.PolynomialRing(
        AbstractAlgebra.QQ, ["x", "y", "z"]
    )
    v = ax > 1 ? ax > 2 ? z : y : x
    po = 0 * v
    if st
        for m in 0:Int(floor(n/2))
            ex = n - 2m
            nm = factorial(big(n)) * (-1)^m
            dn = factorial(big(m)) * factorial(big(ex)) * 2^m
            po += nm//dn * v^ex
        end
    else
        for m in 0:Int(floor(n/2))
            ex = n - 2m
            nm = factorial(big(n)) * (-1)^m
            dn = factorial(big(m)) * factorial(big(ex))
            po += nm//dn * (2v)^ex
        end
    end
    po
end


#------------------------------------------------------------------------------#
#                          Tensor Component Functions                          #
#------------------------------------------------------------------------------#


"""
# Description

    HTC(
        id::String; D::Int = 2, st::Bool = true
    )::AbstractAlgebra.Generic.MPoly{Rational{BigInt}}

Returns a single symbolic *component* of a Hermite Tensor Polynomial.

The indices are specified by `id`, a `String` containing x's, y's,  and/or  z's;
invalid characters on `id` are ignored, while `id`  without  x's,  y's,  or  z's
represent a zero-rank  tensor  (whose  only  component  is  `1`  for  all  space
dimensions).

The Euclidean space dimensionality is specified  by  `D`.  The  tensor  rank  is
implied by the provided `id`, and is always the total  count  of  *valid*  index
characters, considering `D`.

A "probabilist's" polynomial is  returned  if  `st`  is  `true`  (the  default),
otherwise, a "physicist's".

# Usage

This function is not used directly; rather through the interfaces defined on the
`HPT` submodule, which is loaded automatically by the package, by specifying the
AbstractAlgebra ("aa") backend:

```julia-repl
julia> using HermitePolyTensor

julia> HPT.aaHTC.HTC("xxy")
x^2*y-y

```
"""
function HTC(
    id::String; D::Int = 2, st::Bool = true
)::AbstractAlgebra.Generic.MPoly{Rational{BigInt}}
    # Validations
    if D <= 0 || D >= 4
        throw(DomainError("dimension D = $D is outside the valid domain [1, 3]"))
    end
    # Execution
    xc = count(i -> (i == 'x'), id)
    yc = D > 1 ? count(i -> (i == 'y'), id) : 0
    zc = D > 2 ? count(i -> (i == 'z'), id) : 0
    ret = SVHP(xc, ax = 1, st = st)
    if yc != 0
        ret *= SVHP(yc, ax = 2, st = st)
    end
    if zc != 0
        ret *= SVHP(zc, ax = 3, st = st)
    end
    return ret
end


end # module

#------------------------------------------------------------------------------#
#                                    Module                                    #
#------------------------------------------------------------------------------#

"""
# Description

    module seHTC

SymEngine backend (se) module to the Hermite Tensor Polynomial package.

# Usage

This module is `include`d by the `HermitePolyTensor` package and is not  usually
used in isolation, or directly; rather through the  interfaces  defined  on  the
`HPT` submodule, which is loaded automatically by the package:

```julia-repl
julia> using HermitePolyTensor

julia> typeof(HPT.seHTC)
Module

```
"""
module seHTC


#------------------------------------------------------------------------------#
#                                   Imports                                    #
#------------------------------------------------------------------------------#

import SymEngine


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
    )::SymEngine.Basic

Returns a symbolic `n`-th order single-variable Hermite polynomial in  the  axis
`ax`,  axis  1,  2,  and  3  corresponding  to  the  x,  y,  and  z  directions,
respectively. A "probabilist's" polynomial is returned if `st`  is  `true`  (the
default), otherwise, a "physicist's".

This function is not used directly; rather through the interfaces defined on the
`HPT` submodule, which is loaded automatically by the package, by specifying the
SymEngine ("se") backend:

# Usage

```julia-repl
julia> using HermitePolyTensor

julia> HPT.SVHP(3, be = "se")
-3*x + x^3

julia> HPT.SVHP(4, ax = 2, be = "se")
3 - 6*y^2 + y^4

julia> HPT.seHTC.SVHP(5, ax = 3, st = false) # Direct access
120*z - 160*z^3 + 32*z^5

```
"""
function SVHP(
    n::Int; ax::Int = 1, st::Bool = true
)::SymEngine.Basic
    # Validations
    if n < 0
        throw(DomainError("n = $n outside the valid domain [0, ∞)"))
    elseif ax <= 0 || ax >= 4
        throw(DomainError("ax = $ax outside the valid domain [1, 3]"))
    end
    # Execution
    if n == 0
        return SymEngine.Basic(1)
    else
        v = SymEngine.symbols("xyz"[ax:ax])
        ξ = sqrt(v^2) # ξ: U3be
        if st
            ω = (2*SymEngine.PI)^(-1//2) * SymEngine.exp(-ξ^2/2)
        else
            ω = (2*SymEngine.PI)^(-1//2) * SymEngine.exp(-ξ^2)
        end
        r = SymEngine.diff(ω, v, n)
        r *= (-1)^n / ω
        return SymEngine.expand(r)
    end
end


#------------------------------------------------------------------------------#
#                          Tensor Component Functions                          #
#------------------------------------------------------------------------------#


"""
# Description

    HTC(
        id::String; D::Int = 2, st::Bool = true
    )::SymEngine.Basic

Returns a single symbolic *component* of a Hermite tensor polynomial.

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
SymEngine ("se") backend:

```julia-repl
julia> using HermitePolyTensor

julia> HPT.seHTC.HTC("xxy")
-y + x^2*y

```
"""
function HTC(
    id::String; D::Int = 2, st::Bool = true
)::SymEngine.Basic
    # Validations
    if D <= 0 || D >= 4
        throw(DomainError("dimension D = $D outside the valid domain [1, 3]"))
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
    return SymEngine.expand(ret)
end


end # module

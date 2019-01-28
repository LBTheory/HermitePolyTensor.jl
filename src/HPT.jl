#------------------------------------------------------------------------------#
#                                    Module                                    #
#------------------------------------------------------------------------------#

"""
# Description

    module HPT

High-level interface module to the Hermite Polynomial Tensor package.

# Usage

This module is automatically `include`d by the `HermitePolyTensor` package,  and
provides a common interface to the various supported backends.

```julia-repl
julia> using HermitePolyTensor

julia> typeof(HPT)
Module

```
"""
module HPT


#------------------------------------------------------------------------------#
#                                   Imports                                    #
#------------------------------------------------------------------------------#

import Combinatorics


#------------------------------------------------------------------------------#
#                                   Includes                                   #
#------------------------------------------------------------------------------#

include("backend/se/seHTC.jl")
import .seHTC   # Never do "using .seHTC" to avoid function name clashing

include("backend/aa/aaHTC.jl")
import .aaHTC   # Never do "using .aaHTC" to avoid function name clashing


#------------------------------------------------------------------------------#
#                                   Exports                                    #
#------------------------------------------------------------------------------#

# NOTE: This module is @reexport'ed by the main package module

export SVHP, HTC

export HT


#------------------------------------------------------------------------------#
#                                  Constants                                   #
#------------------------------------------------------------------------------#

# Valid Backends
vBE = Dict(
    "aa" => ("AA", "AbsAlg", "AbstractAlgebra"),
    "se" => ("SE", "SymEng", "SymEngine"),
)

# Default Backend
dBE = "aa"  # Since it's way faster than "se"


#------------------------------------------------------------------------------#
#                             Auxiliary Functions                              #
#------------------------------------------------------------------------------#

"""
    theBE(BE::String)::String

Returns the backend name corresponding to `BE`---the backend name or  alias,  or
the default backend `String` for invalid name/alias.
"""
function theBE(BE::String)::String
    for (k, v) in vBE
        if BE == k
            return k
        elseif BE in v
            return k
        end
    end
    return dBE
end


#------------------------------------------------------------------------------#
#                             Interface Functions                              #
#------------------------------------------------------------------------------#

"""
# Description

    SVHP(
        n::Int; ax::Int = 1, st::Bool = true, be::String = dBE
    )

Using  the  backend  specified  by  `be`,  returns  a  symbolic   `n`-th   order
Single-Variable Hermite Polynomial (SVHP) in the axis `ax`, with axis 1, 2,  and
3 corresponding to the x, y, and z directions, respectively.  A  "probabilist's"
polynomial  is  returned  if  `st`  is  `true`  (the  default),   otherwise,   a
"physicist's".

# Usage

```julia-repl
julia> using HermitePolyTensor

julia> SVHP(3, be = "SE")   # SymEngine backend
-3*x + x^3

julia> SVHP(4, ax = 2)      # Default AbstractAlgebra backend
y^4-6//1*y^2+3//1

```
"""
function SVHP(
    n::Int; ax::Int = 1, st::Bool = true, be::String = dBE
)
    # Execution
    BE = theBE(be)
    if BE == "aa"
        return aaHTC.SVHP(n, ax = ax, st = st)
    elseif BE == "se"
        return seHTC.SVHP(n, ax = ax, st = st)
    end
end

"""
# Description

    HTC(
        id::String; D::Int = 2, st::Bool = true, be::String = dBE
    )

Using the backend specified by `be`, returns a symbolic *component* of a Hermite
Tensor Polynomial.

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

```julia-repl
julia> using HermitePolyTensor

julia> HTC("xzyxy", 3)      # Default AbstractAlgebra "aa" backend
x^2*y^2*z-x^2*z-y^2*z+z

julia> HTC("xzyxy", 3, be = "se")   # SymEngine "se" backend
z - x^2*z - y^2*z + x^2*y^2*z

julia> HTC("xzyxy", 3, st = false)  # Physicists's H(*x*)
32//1*x^2*y^2*z-16//1*x^2*z-16//1*y^2*z+8//1*z

```
"""
function HTC(
    id::String; D::Int = 2, st::Bool = true, be::String = dBE
)
    # Execution
    BE = theBE(be)
    if BE == "aa"
        return aaHTC.HTC(id, D = D, st = st)
    elseif BE == "se"
        return seHTC.HTC(id, D = D, st = st)
    end
end


#------------------------------------------------------------------------------#
#                          Multi-Component Functions                           #
#------------------------------------------------------------------------------#

"""
# Description

    HT(
        n::Int; D::Int = 2, st::Bool = true, be::String = dBE
    )::Dict{String,NamedTuple{(:val, :set),Tuple{Any,Set{String}}}}

Using the backend specified by the `be` kwarg, returns  a  dictionary  with  all
*unique* components of a Hermite Tensor Polynomial of order (and rank) `n`, in a
`D`-dimensional Euclidean space.

The dictionary keys are  the  "normalized"  unique  components  of  the  Hermite
Polynomial Tensor, defined by indices `"x"^ax * "y"^ay * "z"^az`,  according  to
the `HTC` function  index  conventions,  that  are  solution  to  the  following
Diophantine [1, 2] equation: `ax + ay + az = n`—there are  `((D  n))`,  or,  `D`
multichoose `n` such indices, or equivalently, `(D + n - 1) choose D` = `(D +  n
- 1)! / (D! * (n - 1)!)` indices alike [2].

Dictionary entries are `namedTuple`'s with symbols (i) `:val` and  (ii)  `:set`,
whose corresponding values are (i) the polynomial component, and (ii) the  index
`Set` to the components having the same  polynomial,  including  the  dictionary
key. Since such `Set` elements are concatenated permutations of  the  normalized
index multiset, with multiplicities `ax`, `ay`, and `az`, the number of multiset
permutations is given by the multinomial coefficient [3] and is equal to  `n!  /
(ax! * ay! * az!)`, or, conversely, by `(ax + ay + az)! / (ax! * ay! * az!)`.

"Probabilist's" polynomials are  returned  if  `st`  is  `true`  (the  default),
otherwise, "physicist's" ones.

# Usage

```julia-repl
```

# References

[1]:  Wikipedia  contributors.  "Diophantine  equation."  Wikipedia,  The   Free
Encyclopedia. Wikipedia, The Free Encyclopedia, 19 Dec. 2018. Web. 25 Jan. 2019.

[2]: Wikipedia contributors. "Combination." Wikipedia,  The  Free  Encyclopedia.
Wikipedia, The Free Encyclopedia, 22 Jan. 2019. Web. 25 Jan. 2019.

[3]: Wikipedia contributors. "Permutation." Wikipedia,  The  Free  Encyclopedia.
Wikipedia, The Free Encyclopedia, 23 Jan. 2019. Web. 25 Jan. 2019.
"""
function HT(
    n::Int; D::Int = 2, st::Bool = true, be::String = dBE
)::Dict{String,NamedTuple{(:val, :set),Tuple{Any,Set{String}}}}
    # Validations
    if n < 0
        throw(DomainError("n = $n outside the valid domain [0, ∞)"))
    elseif D <= 0 || D >= 4
        throw(DomainError("dimension D = $D is outside the valid domain [1, 3]"))
    end
    # Execution
    BE = theBE(be)
    ret = Dict{String,NamedTuple{(:val, :set),Tuple{Any,Set{String}}}}()
    pat = D == 1 ? "x"^n : D == 2 ? "x"^n * "y"^n : "x"^n * "y"^n * "z"^n
    for key in [
        join(i) for i in
        Combinatorics.multiset_combinations(pat, n)
    ]
        theVal = HTC(key, D = D, st = st, be = BE)
        theSet = Set([join(i) for i in
            Combinatorics.multiset_permutations(key, n)
        ])
        ret[key] = (val = theVal, set = theSet)
    end
    ret
end


end # module

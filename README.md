# HermitePolyTensor.jl

Hermite Polynomial Tensor package — Symbolic Hermite Polynomial Tensors

## 1. Description

This package allows for obtaining Hermite Polynomial Tensors — or  N-dimensional
Hermite polynomials [1] — in symbolic form, in  one  to  three  Euclidean  space
dimensions. Currently, it uses two backends, namely, `SymEngine` — identified by
`"se"` — and `AbstractAlgebra` — identified by `"aa"`.

Polynomial types are backend-dependent; they are simply tensor *components*, and
carry no tensor information, nor are able to perform any tensor operation.

Tensor *information* (in a compact form) can nonetheless be  obtained  from  the
package, as a `Dict` type, that behaves as a normal Julia `Dict`,  and  thus  is
unable to perform any tensor operation, but encodes all the  information  needed
to reconstruct the tensor (using external packages).

For   the   `SymEngine`   backend,   tensor   polynomial   components   are   of
`::SymEngine.Basic` type, and `SymEngine` methods  should  be  used  to  further
manipulate them.

For    the    `AbstractAlgebra`    backend,    tensor    components    are    of
`AbstractAlgebra.Generic.MPoly{Rational{BigInt}}`  type,  and  `AbstractAlgebra`
methods should be used to further manipulate the polynomial components.

Of particular interest are  the  so-called  "probabilist's"  Hermite  Polynomial
Tensors [2] — as opposed  to  the  so-called  "physicist's"  Hermite  Polynomial
Tensors [2]. Either can be  obtained,  but  the  "probabilist's"  ones  are  the
default setting for all function methods.

### 1.1. Feature Summary

In summary, with this package, one *can*:

+ Obtain `n`-th order single-variable Hermite polynomials  [2]  symbolically  in
`x`, `y`, or `z`, through the  registered  backends,  with  the  `SVHP`  (Single
Variable Hermite Polynomial) function;

+ Obtain Hermite Polynomial Tensor [1] *components* (one per call)  symbolically
for Euclidean  spaces  of  one  to  three  dimensions,  through  the  registered
backends,  with  the  `HTC`  (Hermite  Tensor  Component)  function,  by   fully
specifying the component index;

+ Obtain a dictionary keyed by tensor indices, each entry containing a  *unique*
(for the given tensor) N-dimensional Hermite Polynomial component, and  a  `Set`
of indices that share the  same  polynomial,  with  the  `HT`  (Hermite  Tensor)
function, by specifying the tensor rank.

In summary, this package (this is not an exaustive list):

- *Does not* implement any tensor container types;

- *Does not* implement any tensor operation functionality;

- *Does not* provide support for Hermite Polynomial Tensor components for spaces
of Euclidean dimension higher than three;

- Although there are tests to the package, it *does not* provide any warranty of
any kind (see the LICENSE file for more details).

## 2. Tutorial

To import the package,  one  can  use  Julia's  standard  `using`  and  `import`
statements. Care has been taken as not to flood the namespace with  symbols,  so
that only submodules and user functions are exported  to  the  target  namespace
(according to Julia's import/using rules).

Exported user functions are: `SVHP`, `HTC`, and `HC` *only*.

```julia-repl
julia> using HermitePolyTensor

julia> HermitePolyTensor.version()  # `version()` is *not* exported
v"0.0.0"

```

The following examples assume a previous `using HermitePolyTensor` statement has
been successfully executed.

### 2.1. Obtaining an `n`-th order single-variable Hermite polynomial

Obtaining a 6-th order, "probabilist" Hermite polynomial  in  `x`  (the  default
axis), i.e., `He₆(x) = x⁶ - 15x⁴ + 45x² - 15` [2]  using  the  `AbstractAlgebra`
backend (the default):

```julia-repl
julia> SVHP(6)
x^6-15//1*x^4+45//1*x^2-15//1

```

Obtaining the previous polynomial using the `SymEngine` backend:

```julia-repl
julia> SVHP(6, be = "se")   # A few aliases for "se" are allowed
-15 + 45*x^2 - 15*x^4 + x^6

julia> HPT.vBE      # valid backends dictionary: name => (alias, alias, ...)
Dict{String,Tuple{String,String,String}} with 2 entries:
  "se" => ("SE", "SymEng", "SymEngine")
  "aa" => ("AA", "AbsAlg", "AbstractAlgebra")

```

Obtaining a 10-th order, "physicist" Hermite polynomial in `z` (the third axis),
i.e., `H₁₀(x) = 1024x¹⁰ - 23040x⁸ + 161280x⁶ - 403200x⁴ + 302400x² - 30240`  [2]
using the `SymEngine` backend:

```julia-repl
julia> SVHP(10, be = "se", ax = 3, st = false)
-30240 + 302400*z^2 - 403200*z^4 + 161280*z^6 - 23040*z^8 + 1024*z^10

```

There is no hardcoded upper limit for  the  Hermite  polynomial  order  in  this
package.

### 2.2. Obtaining Hermite Polynomial Tensor components

Valid indices for the `HTC` function are a concatenation of  `"x"`,  `"y"`,  and
`"z"`, as a Julia `String`, and correspond to the mathematical indices `1`, `2`,
and `3`, respectively, so that `"xxyz"` mean a rank-4 tensor index  `1123`,  for
instance.

Let `H⁽ⁿ⁾(X)`  denote  the  `n`-th  order  (probabilist's)  Hermite  Tensor  for
polynomials in `X = (x[, y[, z]])` defined in a `D`-dimensional Euclidean space,
`1 ≤ D ≤ 3`, and `Hₖₗₘₙ...(X)` (`n` indices `i`  with  `1  ≤  i  ≤  D`)  be  the
`klmn...` polynomial component of the `n`-th order Hermite polynomial  tensor  —
in general, a multivariate polynomial in `X`.

Obtaining  the  `H₁₂₁₁₂₁₂(X)`  polynomial  component  of  `H⁽⁷⁾(X)`  for  a  2-D
Euclidean space, i.e., the `-9y +18x²y -6x²y³ -3x⁴y³ +3y³` polynomial, using the
default `AbstractAlgebra` backend, and then the `SymEngine` backend, is done  as
follows:

```julia-repl
julia> HTC("xyxxyxy", 2)
x^4*y^3-3//1*x^4*y-6//1*x^2*y^3+18//1*x^2*y+3//1*y^3-9//1*y

julia> HTC("xyxxyxy", 2, be = "SymEngine")
-9*y + 18*x^2*y - 6*x^2*y^3 - 3*x^4*y + x^4*y^3 + 3*y^3

```

One property of Hermite Polynomial Tensors is that all components  described  by
permutations of a given index set or multiset are the same, so that  `H₁₂₁(X)  =
H₁₁₂(X) = H₂₁₁(X)`, for instance:

```julia-repl
julia> HTC("xyx") == HTC("xxy") == HTC("yxx")
true

```

Therefore, high-rank Hermite Polynomial Tensor components can  be  specified  by
string repetition and concatenation:

```julia-repl
julia> SVHP(9, be = "se")
945*x - 1260*x^3 + 378*x^5 - 36*x^7 + x^9

julia> HTC("x"^9, be = "se")
945*x - 1260*x^3 + 378*x^5 - 36*x^7 + x^9

julia> HTC("x"^6 * "y"^3, be = "se")
45*y - 135*x^2*y + 45*x^2*y^3 + 45*x^4*y - 15*x^4*y^3 - 3*x^6*y + x^6*y^3 - 15*y^3

```

### 2.3. Obtaining Hermite Polynomial Tensor information

The `HT` function returns a dictionary that encodes  Hermite  Polynomial  Tensor
data, avoiding repeats:

Obtaining the 3rd-order Hermite Polynomial Tensor information in  two  Euclidean
space dimensions:

```julia-repl
julia> HT(3)
Dict{String,NamedTuple{(:val, :set),Tuple{Any,Set{String}}}} with 4 entries:
  "xxx" => NamedTuple{(:val, :set),Tuple{Any,Set{String}}}((x^3-3//1*x, Set(["xxx"])))
  "xyy" => NamedTuple{(:val, :set),Tuple{Any,Set{String}}}((x*y^2-x, Set(["xyy", "yyx", "yxy"])))
  "yyy" => NamedTuple{(:val, :set),Tuple{Any,Set{String}}}((y^3-3//1*y, Set(["yyy"])))
  "xxy" => NamedTuple{(:val, :set),Tuple{Any,Set{String}}}((x^2*y-y, Set(["yxx", "xyx", "xxy"])))

```

The output shows that the `"xyy"` component has a multiplicity of  three  —  the
`x*y^2-x` polynomial is found on components `"xyy"`, `"yyx"`, and  `"yxy"`  (all
members of the `Set`).

Obtaining the 4th-order Hermite Polynomial Tensor information in three Euclidean
space dimensions, using the `"SymEngine"` (`"se"`) backend:

```julia-repl
julia> HT(4, D = 3, be = "se")
Dict{String,NamedTuple{(:val, :set),Tuple{Any,Set{String}}}} with 15 entries:
  "yzzz" => NamedTuple{(:val, :set),Tuple{Any,Set{String}}}((-3*y*z + y*z^3, Set(["yzzz", "zzyz", "zzzy", "zyzz"])))
  "xxyz" => NamedTuple{(:val, :set),Tuple{Any,Set{String}}}((-y*z + x^2*y*z, Set(["zxxy", "yxxz", "xxzy", "xxyz", "xyxz", "yzxx…
  "xxzz" => NamedTuple{(:val, :set),Tuple{Any,Set{String}}}((1 + x^2*z^2 - x^2 - z^2, Set(["zzxx", "xzzx", "xxzz", "xzxz", "zxz…
  "yyyy" => NamedTuple{(:val, :set),Tuple{Any,Set{String}}}((3 - 6*y^2 + y^4, Set(["yyyy"])))
  "xxxx" => NamedTuple{(:val, :set),Tuple{Any,Set{String}}}((3 - 6*x^2 + x^4, Set(["xxxx"])))
  "yyyz" => NamedTuple{(:val, :set),Tuple{Any,Set{String}}}((-3*y*z + y^3*z, Set(["yyyz", "zyyy", "yzyy", "yyzy"])))
  "xxxz" => NamedTuple{(:val, :set),Tuple{Any,Set{String}}}((-3*x*z + x^3*z, Set(["zxxx", "xxxz", "xzxx", "xxzx"])))
  "xzzz" => NamedTuple{(:val, :set),Tuple{Any,Set{String}}}((-3*x*z + x*z^3, Set(["zzxz", "zxzz", "xzzz", "zzzx"])))
  "xyyz" => NamedTuple{(:val, :set),Tuple{Any,Set{String}}}((-x*z + x*y^2*z, Set(["yxyz", "zyxy", "xyzy", "zyyx", "xyyz", "yxzy…
  "yyzz" => NamedTuple{(:val, :set),Tuple{Any,Set{String}}}((1 + y^2*z^2 - y^2 - z^2, Set(["yyzz", "zzyy", "yzyz", "yzzy", "zyy…
  "xyyy" => NamedTuple{(:val, :set),Tuple{Any,Set{String}}}((-3*x*y + x*y^3, Set(["xyyy", "yxyy", "yyxy", "yyyx"])))
  "xxyy" => NamedTuple{(:val, :set),Tuple{Any,Set{String}}}((1 + x^2*y^2 - x^2 - y^2, Set(["xyxy", "yxxy", "xxyy", "xyyx", "yyx…
  "xxxy" => NamedTuple{(:val, :set),Tuple{Any,Set{String}}}((-3*x*y + x^3*y, Set(["xyxx", "xxxy", "xxyx", "yxxx"])))
  "zzzz" => NamedTuple{(:val, :set),Tuple{Any,Set{String}}}((3 - 6*z^2 + z^4, Set(["zzzz"])))
  "xyzz" => NamedTuple{(:val, :set),Tuple{Any,Set{String}}}((-x*y + x*y*z^2, Set(["zyzx", "zzxy", "yzxz", "zyxz", "xzyz", "zxzy…

```

The output shows that  the  tensor  is  comprised  of  15  different  polynomial
components (15 entries to the `Dict` output). All such polynomials are displayed
(see more below on how to get a specific component).

Probing the previous output:

```julia-repl
julia> HT(4, D = 3, be = "se")["xxzz"].val
1 + x^2*z^2 - x^2 - z^2

julia> HT(4, D = 3, be = "se")["xxzz"].set
Set(["zzxx", "xzzx", "xxzz", "xzxz", "zxzx", "zxxz"])

```

## 3. Benchmarks

Low-rank Hermite polynomial tensor component by backend:

```julia-repl
julia> using BenchmarkTools

julia> @benchmark HTC("xxy", be = "aa")
BenchmarkTools.Trial:
  memory estimate:  21.78 KiB
  allocs estimate:  765
  --------------
  minimum time:     28.571 μs (0.00% GC)
  median time:      32.596 μs (0.00% GC)
  mean time:        70.327 μs (39.18% GC)
  maximum time:     111.622 ms (81.25% GC)
  --------------
  samples:          10000
  evals/sample:     1

julia> @benchmark HTC("xxy", be = "se")
BenchmarkTools.Trial:
  memory estimate:  11.62 KiB
  allocs estimate:  1274
  --------------
  minimum time:     98.678 μs (0.00% GC)
  median time:      101.177 μs (0.00% GC)
  mean time:        129.467 μs (6.80% GC)
  maximum time:     81.879 ms (60.12% GC)
  --------------
  samples:          10000
  evals/sample:     1

```

High-rank Hermite polynomial tensor component by backend:

```julia-repl
julia> @benchmark HTC("x"^20 * "y"^20 * "z"^20, be = "aa")
BenchmarkTools.Trial:
  memory estimate:  207.45 KiB
  allocs estimate:  8943
  --------------
  minimum time:     334.628 μs (0.00% GC)
  median time:      352.372 μs (0.00% GC)
  mean time:        750.793 μs (33.74% GC)
  maximum time:     134.450 ms (77.26% GC)
  --------------
  samples:          6640
  evals/sample:     1

julia> @benchmark HTC("x"^20 * "y"^20 * "z"^20, be = "se")
BenchmarkTools.Trial:
  memory estimate:  323.69 KiB
  allocs estimate:  41154
  --------------
  minimum time:     3.088 ms (0.00% GC)
  median time:      3.116 ms (0.00% GC)
  mean time:        3.495 ms (1.50% GC)
  maximum time:     69.311 ms (64.29% GC)
  --------------
  samples:          1429
  evals/sample:     1

```

## References

[1] H. Grad, “Note on N-dimensional Hermite polynomials,” Communications on Pure
and Applied Mathematics, vol. 2, no. 4, pp. 325–330, 1949.

[2]  Wikipedia  contributors.  “Hermite  polynomials.”   Wikipedia,   The   Free
Encyclopedia. Wikipedia, The Free Encyclopedia, 9 Jan. 2019. Web. 23 Jan. 2019.

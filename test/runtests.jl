#------------------------------------------------------------------------------#
#                                   Includes                                   #
#------------------------------------------------------------------------------#

using Test
import HermitePolyTensor
import SymEngine
import Random
import AbstractAlgebra

#------------------------------------------------------------------------------#
#     Hermite Tensor Polynomial Component (HTC) Tests -- SymEngine backend     #
#------------------------------------------------------------------------------#

# Preparations
# ------------

x, y, z = SymEngine.symbols("x y z")

He3x = x^3 - 3x
He6x = x^6 - 15x^4 + 45x^2 - 15
He9x = x^9 - 36x^7 + 378x^5 - 1260x^3 + 945x

He3y = y^3 - 3y
He6y = y^6 - 15y^4 + 45y^2 - 15
He9y = y^9 - 36y^7 + 378y^5 - 1260y^3 + 945y

He3z = z^3 - 3z
He6z = z^6 - 15z^4 + 45z^2 - 15
He9z = z^9 - 36z^7 + 378z^5 - 1260z^3 + 945z

H_3x = 8x^3 - 12x
H_6x = 64x^6 - 480x^4 + 720x^2 - 120
H_9x = 512x^9 - 9216x^7 + 48384x^5 - 80640x^3 + 30240x

H_3y = 8y^3 - 12y
H_6y = 64y^6 - 480y^4 + 720y^2 - 120
H_9y = 512y^9 - 9216y^7 + 48384y^5 - 80640y^3 + 30240y

H_3z = 8z^3 - 12z
H_6z = 64z^6 - 480z^4 + 720z^2 - 120
H_9z = 512z^9 - 9216z^7 + 48384z^5 - 80640z^3 + 30240z


# Single-variable HTC Tests
# -------------------------

@testset "SVHP, SymEng, x-axis, prob. Hermite poly          " begin
    @test He3x == HermitePolyTensor.SVHP(3, be = "se", ax = 1, st = true)
    @test He6x == HermitePolyTensor.SVHP(6, be = "se", ax = 1, st = true)
    @test He9x == HermitePolyTensor.SVHP(9, be = "se", ax = 1, st = true)
end

@testset "SVHP, SymEng, y-axis, prob. Hermite poly          " begin
    @test He3y == HermitePolyTensor.SVHP(3, be = "se", ax = 2, st = true)
    @test He6y == HermitePolyTensor.SVHP(6, be = "se", ax = 2, st = true)
    @test He9y == HermitePolyTensor.SVHP(9, be = "se", ax = 2, st = true)
end

@testset "SVHP, SymEng, z-axis, prob. Hermite poly          " begin
    @test He3z == HermitePolyTensor.SVHP(3, be = "se", ax = 3, st = true)
    @test He6z == HermitePolyTensor.SVHP(6, be = "se", ax = 3, st = true)
    @test He9z == HermitePolyTensor.SVHP(9, be = "se", ax = 3, st = true)
end

@testset "SVHP, SymEng, x-axis, phys. Hermite poly          " begin
    @test H_3x == HermitePolyTensor.SVHP(3, be = "se", ax = 1, st = false)
    @test H_6x == HermitePolyTensor.SVHP(6, be = "se", ax = 1, st = false)
    @test H_9x == HermitePolyTensor.SVHP(9, be = "se", ax = 1, st = false)
end

@testset "SVHP, SymEng, y-axis, phys. Hermite poly          " begin
    @test H_3y == HermitePolyTensor.SVHP(3, be = "se", ax = 2, st = false)
    @test H_6y == HermitePolyTensor.SVHP(6, be = "se", ax = 2, st = false)
    @test H_9y == HermitePolyTensor.SVHP(9, be = "se", ax = 2, st = false)
end

@testset "SVHP, SymEng, z-axis, phys. Hermite poly          " begin
    @test H_3z == HermitePolyTensor.SVHP(3, be = "se", ax = 3, st = false)
    @test H_6z == HermitePolyTensor.SVHP(6, be = "se", ax = 3, st = false)
    @test H_9z == HermitePolyTensor.SVHP(9, be = "se", ax = 3, st = false)
end

# Hermite Tensor Polynomial Component Tests
#------------------------------------------

@testset "HTC, SymEng, Factorial, 3D, prob. Hermite T. Comp." begin
    for X in [(He3x, 3), (He6x, 6), (He9x, 9)]
        for Y in [(He3y, 3), (He6y, 6), (He9y, 9)]
            for Z in [(He3z, 3), (He6z, 6), (He9z, 9)]
                raw = "x"^X[2] * "y"^Y[2] * "z"^Z[2]    # raw (ordered) index set
                idx = join(Random.shuffle([i for i in raw]))   # shuffled index set
                @test SymEngine.expand(X[1] * Y[1] * Z[1]) == HermitePolyTensor.HTC(idx, D = 3, be = "se", st = true)
            end
        end
    end
end

@testset "HTC, SymEng, Factorial, 3D, phys. Hermite T. Comp." begin
    for X in [(H_3x, 3), (H_6x, 6), (H_9x, 9)]
        for Y in [(H_3y, 3), (H_6y, 6), (H_9y, 9)]
            for Z in [(H_3z, 3), (H_6z, 6), (H_9z, 9)]
                raw = "x"^X[2] * "y"^Y[2] * "z"^Z[2]    # raw (ordered) index set
                idx = join(Random.shuffle([i for i in raw]))   # shuffled index set
                @test SymEngine.expand(X[1] * Y[1] * Z[1]) == HermitePolyTensor.HTC(idx, D = 3, be = "se", st = false)
            end
        end
    end
end


#------------------------------------------------------------------------------#
#  Hermite Tensor Polynomial Component (HTC) Tests -- AbstractAlgebra backend  #
#------------------------------------------------------------------------------#

# Preparations
# ------------

R, (x, y, z) = AbstractAlgebra.PolynomialRing(AbstractAlgebra.QQ, ["x", "y", "z"])

He3x = x^3 - 3x
He6x = x^6 - 15x^4 + 45x^2 - 15
He9x = x^9 - 36x^7 + 378x^5 - 1260x^3 + 945x

He3y = y^3 - 3y
He6y = y^6 - 15y^4 + 45y^2 - 15
He9y = y^9 - 36y^7 + 378y^5 - 1260y^3 + 945y

He3z = z^3 - 3z
He6z = z^6 - 15z^4 + 45z^2 - 15
He9z = z^9 - 36z^7 + 378z^5 - 1260z^3 + 945z

H_3x = 8x^3 - 12x
H_6x = 64x^6 - 480x^4 + 720x^2 - 120
H_9x = 512x^9 - 9216x^7 + 48384x^5 - 80640x^3 + 30240x

H_3y = 8y^3 - 12y
H_6y = 64y^6 - 480y^4 + 720y^2 - 120
H_9y = 512y^9 - 9216y^7 + 48384y^5 - 80640y^3 + 30240y

H_3z = 8z^3 - 12z
H_6z = 64z^6 - 480z^4 + 720z^2 - 120
H_9z = 512z^9 - 9216z^7 + 48384z^5 - 80640z^3 + 30240z


# Single-variable HTC Tests
# -------------------------

@testset "SVHP, AbsAlg, x-axis, prob. Hermite poly          " begin
    @test He3x == HermitePolyTensor.SVHP(3, be = "aa", ax = 1, st = true)
    @test He6x == HermitePolyTensor.SVHP(6, be = "aa", ax = 1, st = true)
    @test He9x == HermitePolyTensor.SVHP(9, be = "aa", ax = 1, st = true)
end

@testset "SVHP, AbsAlg, y-axis, prob. Hermite poly          " begin
    @test He3y == HermitePolyTensor.SVHP(3, be = "aa", ax = 2, st = true)
    @test He6y == HermitePolyTensor.SVHP(6, be = "aa", ax = 2, st = true)
    @test He9y == HermitePolyTensor.SVHP(9, be = "aa", ax = 2, st = true)
end

@testset "SVHP, AbsAlg, z-axis, prob. Hermite poly          " begin
    @test He3z == HermitePolyTensor.SVHP(3, be = "aa", ax = 3, st = true)
    @test He6z == HermitePolyTensor.SVHP(6, be = "aa", ax = 3, st = true)
    @test He9z == HermitePolyTensor.SVHP(9, be = "aa", ax = 3, st = true)
end

@testset "SVHP, AbsAlg, x-axis, phys. Hermite poly          " begin
    @test H_3x == HermitePolyTensor.SVHP(3, be = "aa", ax = 1, st = false)
    @test H_6x == HermitePolyTensor.SVHP(6, be = "aa", ax = 1, st = false)
    @test H_9x == HermitePolyTensor.SVHP(9, be = "aa", ax = 1, st = false)
end

@testset "SVHP, AbsAlg, y-axis, phys. Hermite poly          " begin
    @test H_3y == HermitePolyTensor.SVHP(3, be = "aa", ax = 2, st = false)
    @test H_6y == HermitePolyTensor.SVHP(6, be = "aa", ax = 2, st = false)
    @test H_9y == HermitePolyTensor.SVHP(9, be = "aa", ax = 2, st = false)
end

@testset "SVHP, AbsAlg, z-axis, phys. Hermite poly          " begin
    @test H_3z == HermitePolyTensor.SVHP(3, be = "aa", ax = 3, st = false)
    @test H_6z == HermitePolyTensor.SVHP(6, be = "aa", ax = 3, st = false)
    @test H_9z == HermitePolyTensor.SVHP(9, be = "aa", ax = 3, st = false)
end

# Hermite Tensor Polynomial Component Tests
#------------------------------------------

@testset "HTC, AbsAlg, Factorial, 3D, prob. Hermite T. Comp." begin
    for X in [(He3x, 3), (He6x, 6), (He9x, 9)]
        for Y in [(He3y, 3), (He6y, 6), (He9y, 9)]
            for Z in [(He3z, 3), (He6z, 6), (He9z, 9)]
                raw = "x"^X[2] * "y"^Y[2] * "z"^Z[2]    # raw (ordered) index set
                idx = join(Random.shuffle([i for i in raw]))   # shuffled index set
                @test (X[1] * Y[1] * Z[1]) == HermitePolyTensor.HTC(idx, D = 3, be = "aa", st = true)
            end
        end
    end
end

@testset "HTC, AbsAlg, Factorial, 3D, phys. Hermite T. Comp." begin
    for X in [(H_3x, 3), (H_6x, 6), (H_9x, 9)]
        for Y in [(H_3y, 3), (H_6y, 6), (H_9y, 9)]
            for Z in [(H_3z, 3), (H_6z, 6), (H_9z, 9)]
                raw = "x"^X[2] * "y"^Y[2] * "z"^Z[2]    # raw (ordered) index set
                idx = join(Random.shuffle([i for i in raw]))   # shuffled index set
                @test (X[1] * Y[1] * Z[1]) == HermitePolyTensor.HTC(idx, D = 3, be = "aa", st = false)
            end
        end
    end
end


#------------------------------------------------------------------------------#
#                     Hermite Polynomial Tensor (HT) Tests                     #
#------------------------------------------------------------------------------#

function nChooseK(n, k)
    return Int(factorial(big(n)) / (factorial(big(k)) * factorial(big(n - k))))
end

function nMultiChooseK(n, k)
    return nChooseK(n + k - 1, n)
end

@testset "HT, unique tensor components                      " begin
    for theD in 1:3
        for theN in 1:10
            @test length(HermitePolyTensor.HT(theN, D = theD)) == nMultiChooseK(theN, theD)
        end
    end
end

function nMultiPerms(M)
    return Int(factorial(big(sum(M))) / prod(map(factorial, map(big, M))))
end

@testset "HT, unique component multiplicities               " begin
    for theD in 1:3
        for theN in 1:10
            H = HermitePolyTensor.HT(theN, D = theD)
            for id in keys(H)
                xc = count(i -> (i == 'x'), id)
                yc = theD > 1 ? count(i -> (i == 'y'), id) : 0
                zc = theD > 1 ? count(i -> (i == 'z'), id) : 0
                @test length(H[id].set) == nMultiPerms((xc, yc, zc))
            end
        end
    end
end

@testset "HT, total tensor components                       " begin
    for theD in 1:3
        for theN in 1:10
            @test sum(length(k.set) for k in values(HermitePolyTensor.HT(theN, D = theD))) == theD^theN
        end
    end
end



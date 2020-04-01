module ReflecBasis
    ####################################
        function FriedrichsAngleAB(A::Matrix{Float64},B::Matrix{Float64})
            # Calculating the Friederich Angle
            QA, RA = qr(A')
            QB, RB = qr(B')
            S =  svd(QA'*QB)
            # Angle in Degrees
            # angleAB = acosd(maximum(S))
            # Angle in Radians
            # angleAB = acos(maximum(S))
            # println(S[2])
            ind = findfirst(x -> x<(1-1e-8),S[2])
            # return maximum(S[2])
            return acos(S[2][ind])
        end
    ####################################
    ####################################
        function contructProjector2(A::Matrix{Float64},a::Vector{Float64},n::Int64)
            MA = pinv(A)
            PA=eye(n) - MA*A
            aP = MA*a
            return PA, aP
        end
    ####################################
        function contructProjector(A::Matrix{Float64},a::Vector{Float64},n::Int64)
            QA, RA = qr(A')
            PA= eye(n) - QA*QA'
            aP = QA*(RA'\a)
            return PA, aP
        end
    ####################################
        # function contructProjector(A::SparseMatrixCSC{Float64,Int64},a::Vector{Float64})
        #     MA = pinv(A)
        #     PA=eye(n) - MA*A
        #     aP = MA*a
        #     return PA, aP
        # end
    ####################################
    ####################################
        function Projection(PA::Matrix{Float64},aP::Vector{Float64},xzero::Vector{Float64})
            # Orthogonal projection into affine subspace Ax=b if m<n.
            return PA*xzero + aP

        end
    ####################################
    ####################################
        function Reflection(PA::Matrix{Float64},aP::Vector{Float64},xzero::Vector{Float64})
            # Reflection from affine subespace Ax=b if m<n.
            # using LinearLeastSquares: min||x-xzero|| s.t. Ax=b.
            # Rbar(xzero) = R(xzero-b) + b
            # R = 2P - I where P is a projection
            proj_xzero = Projection(PA,aP,xzero);
            return  2*proj_xzero - xzero

        end
    ####################################

    export FriedrichsAngleAB, GenerateRdmPair, contructProjector, Projection         

end
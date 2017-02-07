__precompile__()

module MAP
    export GenerateSamples, FriedrichsAngleAB, algorithmMAP, fourpointscheme, algorithmDRM
    export algorithmDRM_C, GenerateRdmPair, friedrichs

    # global const EPS_VAL = 1e-3
    global const ITER_MAX = 5000000

    # println("Tolerance: $EPS_VAL")
####################################
####################################
    function GenerateSamples(n::Int64,affine::Bool=true)
        ma = rand(1:n-1) ## number of extra normals of A
        mb = rand(1:n-1-ma) ## number of extra normals of B
        A = randn(ma,n)
        B = randn(mb,n)
        if affine
            a = randn(ma)
            b = randn(mb)
        else
            a = zeros(ma)
            b = zeros(mb)
        end

        return A, a, ma, B, b, mb
    end
####################################
####################################
    function GenerateRdmPair(n::Int64,cmax::Int64,affine::Bool=true)
        n >= cmax ||  error("n must grater then common normals")
        mcex = rand(1:cmax) ## number of common normals
        maex = rand(1:n-1-mcex) ## number of extra normals of A
        mbex = rand(1:n-1-mcex-maex) ## number of extra normals of B


        Cex = randn(mcex,n)
        while !(rank(Cex) == mcex)
            Cex = randn(mcex,n)
        end

        Aex = randn(maex,n)
        while !(rank([Aex; Cex]) == maex +  mcex)
            Aex = randn(maex,n)
        end

        Bex = randn(mbex,n)
        while !(rank([Aex; Cex; Bex]) == maex +  mcex + mbex)
            Bex = randn(mbex,n)
        end

        ma = maex + mcex
        mb = mbex + mcex

        # Space U
        A = [Aex; Cex]

        # Space V
        B = [Cex; Bex]

        if affine
            a = randn(ma)
            b = randn(mb)
        else
            a = zeros(ma)
            b = zeros(mb)
        end
        return A, a, ma, B, b, mb
    end
####################################

####################################
function friedrichs(A,B)

    cdm = rank(A)+rank(B)-rank([A;B]);
    C=[A;B[cdm+1:rank(B),:]];

    AA=pinv(A)*A;
    BB=pinv(B)*B;
    CC=pinv(C)*C

    M=CC+BB*AA-BB-AA;

    return sqrt(eigmax(M*M'));

end
####################################

####################################
    function FriedrichsAngleAB(A::Matrix{Float64},B::Matrix{Float64})
        # Calculating Friederich Angle
        QA, RA = qr(A')
        QB, RB = qr(B')
        S =  svd(QA'*QB)
        # Angle in Degrees
        # angleAB = acosd(maximum(S))
        # Angle in Radians
        # angleAB = acos(maximum(S))
        println(S[2])
        ind = findfirst(x -> x<(1-1e-8),S[2])
        # return maximum(S[2])
        return S[2][ind]
    end
####################################
####################################
    function contructProjector(A::Matrix{Float64},a::Vector{Float64},n::Int64)
        MA = pinv(A)
        PA=eye(n) - MA*A
        aP = MA*a
        return PA, aP
    end
####################################
####################################
    function Projection(PA::Matrix{Float64},aP::Vector{Float64},xzero::Vector{Float64})
        # Orthogonal projection into affine subespace Ax=b if m<n.
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
        return xzero + 2*(proj_xzero - xzero)

    end
####################################
####################################
    function  fourpointscheme(A::Matrix{Float64},a::Vector{Float64},B::Matrix{Float64},
                                        b::Vector{Float64},n::Int64,xzero::Vector{Float64},
                                        file::IOStream,printfile=true,method::Int64=1, EPS_VAL::Float64=1e-3)
    #algoritmoRBLRS using 4 points in inner interaction
       #begining the projection onto A and onto B
        PA, aP = contructProjector(A,a,n)
        PB, bP = contructProjector(B,b,n)
        xstarA = Projection(PA,aP,xzero)
        xstarB = Projection(PB,bP,xzero)
        tol_xstarA = norm([A; B]*xstarA - [a; b],2)
        tol_xstarB = norm([A; B]*xstarB - [a; b],2)
        if tol_xstarA < tol_xstarB
            xstar = xstarA
            ind = 1
        else
            xstar = xstarB
            ind = 2
        end
    #   use ind = randi([1 2]) to implement random start
        tol = 1
        k =2
        while (tol >= EPS_VAL)  || (k>ITER_MAX)
            if (mod(ind,2) == 0)
                    b1 = xstar
                    #projection onto A
                    a1 = Projection(PA,aP,b1)
                    #projection onto B
                    b2 = Projection(PB,bP,a1)
                    #projection onto A
                    a2 = Projection(PA,aP,b2)

            else
                    a1 = xstar
                    #projection onto B
                    b1 = Projection(PB,bP,a1)
                    #projection onto A
                    a2 = Projection(PA,aP,b1)
                    #projection onto B
                    b2 = Projection(PB,bP,a2)
            end
            if method == 1
                # No Acceleration
                atil = a2
                btil = b2
            elseif method == 2
                # Closest Point Acceleration
                atil = accCP(a1,b1,a2)
                btil = accCP(b1,a2,b2)
            elseif method == 3
                atil = accEMAP(a1,b1,a2)
                btil = accEMAP(b1,a2,b2)
            end

            # Lembrete: Em Omega Geral, projetar atil/btil em OmA ou OmB
            tol_A = norm([A; B]*a2 - [a; b],2)
            tol_B = norm([A; B]*b2 - [a; b],2)
            tolac_A = norm([A; B]*atil - [a; b],2)
            tolac_B = norm([A; B]*btil - [a; b],2)
            vectors = [a2 b2 atil btil]
            tol, ind = findmin([tol_A tol_B tolac_A tolac_B])
            xstar = vectors[:,ind]
            k+=3
            (k>ITER_MAX) && (warn("Maximum number of projections achievied in 4point"); break)
        end
        if printfile
                str = @sprintf("    conv %7d  %10.8e",k,tol)
                print(file,str)
        else
            println("Number of projections: $k")
            @printf("Error = %s\n", tol)
        end
        # println(xstar)
    end
####################################
####################################
    function  algorithmMAP(A::Matrix{Float64},a::Vector{Float64},B::Matrix{Float64},
                                        b::Vector{Float64},n::Int64,xzero::Vector{Float64},
                                        file::IOStream,printfile=true,method::Int64=1, EPS_VAL::Float64=1e-3)
        # Projecting onto A first
        PA, aP = contructProjector(A,a,n)
        PB, bP = contructProjector(B,b,n)
        xstar = Projection(PA,aP,xzero)
        tol = 1
        k =1
        while (tol >= EPS_VAL)
            a1 = xstar
            #projection onto B
            b1 = Projection(PB,bP,a1)
            #projection onto A
            a2 = Projection(PA,aP,b1)
            if method == 1
                # RAS Acceleration
                xstar = a2
            elseif method == 2
                # Closest Point Acceleration
                xstar = accCP(a1,b1,a2)
            elseif method == 3
                xstar = accEMAP(a1,b1,a2)
            end
            tol = norm([A;B]*xstar - [a;b],2)
            k+=2
            (k>ITER_MAX) && (warn("Maximum number of projections achievied in MAP"); break)
        end
        if printfile
                str = @sprintf("    conv %7d  %10.8e",k,tol)
                print(file,str)
        else
            println("Number of projections: $k")
            @printf("Error = %s\n", tol)
        end
    end
####################################
####################################
    function algorithmDRM(A::Matrix{Float64},a::Vector{Float64},B::Matrix{Float64},
                                        b::Vector{Float64},n::Int64,xzero::Vector{Float64},
                                        file::IOStream,printfile=true,EPS_VAL::Float64=1e-3)
        PA, aP = contructProjector(A,a,n)
        PB, bP = contructProjector(B,b,n)
        xstar = xzero
        k = 0
        tol = norm([A;B]*xstar - [a;b],2)
        while (tol >= EPS_VAL)
            reflec = Reflection(PA,aP,xstar)
            reflec = Reflection(PB,bP,reflec)
            xstar = (xstar + reflec)/2;
            tol = norm([A;B]*xstar - [a;b],2)
            k +=2
            (k>ITER_MAX) && (warn("Maximum number of projections achievied in DRM"); break)
        end
        if printfile
                str = @sprintf("    conv %7d  %10.8e",k,tol)
                print(file,str)
        else
            println("Number of projections: $k")
            @printf("Error = %s\n", tol)
        end
    end
####################################
####################################
    function algorithmDRM_C(A::Matrix{Float64},a::Vector{Float64},
                            B::Matrix{Float64},b::Vector{Float64},
                            n::Int64,xzero::Vector{Float64},file::IOStream,
                            printfile=true,EPS_VAL::Float64=1e-3)
        PA, aP = contructProjector(A,a,n)
        PB, bP = contructProjector(B,b,n)
        xstar = xzero
        k = 0
        tol = 1.0
        while (tol >= EPS_VAL)
            ypto = Reflection(PA,aP,xstar)
            zpto = Reflection(PB,bP,ypto)
            if norm(ypto - xstar)<EPS_VAL
                xstar = 0.5(xstar + zpto)
                tol = norm([A;B]*xstar - [a;b],2)
                continue
            elseif norm(zpto - ypto)<EPS_VAL
                xstar = 0.5(xstar + zpto)
                tol = norm([A;B]*xstar - [a;b],2)
                continue
            end
            v_xy = ypto - xstar
            v_yz = zpto - ypto
            med_xy = xstar + 0.5*(v_xy)
            med_yz = ypto + 0.5*(v_yz)
            dir_xy = v_yz - (dot(v_yz,v_xy)/dot(v_xy,v_xy))*v_xy
            dir_yz = v_xy - (dot(v_xy,v_yz)/dot(v_yz,v_yz))*v_yz
            t = [dir_xy -dir_yz]\(med_yz - med_xy)
            xstar = med_xy + t[1]*dir_xy
            tol = norm([A;B]*xstar - [a;b],2)
            k += 2
            println("Number of projections: $k")
            @printf("Error = %s\n", tol)
            (k>ITER_MAX) && (warn("Maximum number of projections achievied in DRM-C"); break)
        end
        if printfile
                str = @sprintf("    conv %7d  %10.8e",k,tol)
                print(file,str)
        else
            println("Number of projections: $k")
            @printf("Error = %s\n", tol)
        end
    end
####################################
####################################
    function accEMAP(x1,x2,x3)
        # x1 and x3 are in subspace A
        # x3 is in subspace B.
        # Finds the Extrapolation VonNewmann-Accelerated
        # point xacc in subspace given by span(x3-x1)
        t = norm(x1-x2)^2/norm(x3-x1)^2
        return x1 + (x3 - x1)*t
    end
####################################
####################################
    function accCP(x1,x2,x3)
    # Finds the Closest Point Acceleration xac
    # x2 = PB*x1 in B if x1 in A.
    # x3 = PA*PB*x1, if x1 in A.

        # if (x1 ~= x3)
            t = dot(x1,x1-x3)/((norm(x1 - x3,2))^2)
        # else
            # t = 1
        # end
        return x1 + t*(x3 -x1)

    end
####################################

####################################

end

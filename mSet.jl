__precompile__()
# module mSet
    include("CRM.jl")
    if VERSION >=  v"0.7.0"
        import SparseArrays
        import LinearAlgebra
        import LinearAlgebra: qr
        import SuiteSparse
        import Printf
        import Random
    end
    # export algorithmCDRMmSet, GenerateSampleApart

    # global const ZERO_VAL = 1e-15
    # global const iter_max = 1e10


    # println("Tolerance: $EPS_VAL")
####################################
####################################
    function GenerateSampleApart(n::Int64,m::Int64)
        APartition = rand(1:n,m,1)
        R = sum(APartition)
        AA = []
        for I = 1:m
            push!(AA,randn(APartition[I],n))
        end
        return AA, R
    end
####################################




####################################
    function algorithmCDRMmSet(n::Int64,m::Int64,PA,aP,xstar::Vector{Float64},
        xbar::Vector{Float64},#\file::IOStream,
        printfile=false,EPS_VAL::Float64=1e-3)
        tol = 1.0
        conv =  String("conv")
        # xstar =  Projection(PB,bP,xstar)  # Initial Point on span(U,V)
        k = 1
        rate = 0.
        # X is the matrix of the reflections
        X = []
        xstar = Projection(PA[2],aP[2],xstar)
        # println(xstar)
        while (tol >= EPS_VAL)
            push!(X,xstar)
            # xbar = xstar #Looking for FixingPoint
            # println("Entrou no loop $k")
            for J = 1:m
                push!(X,MAP.Reflection(PA[J],aP[J],X[J]))
                # println(X[:,J+1])
            end
            # println("Criou X")
            xstar = FindCircumcentermSet(X)
            # println(xstar)
            # gapA = Projection(PA,aP,xstar)
            # gapB = Projection(PB,bP,xstar)
            # tol = norm(gapA - gapB,2) #Gap Distance
            tol = norm(xstar - xbar,2) #True Error
            # println(tol)
            # rate = tol/norm(xold - xbar,2)
            k += m
            (k>iter_max) && (warn("Maximum number of projections achieved in C-DRM"); conv=String("itmax"); break)
        end
        # println(rate)
        if printfile
                str = @sprintf("    %s  %7d  %10.8e %10.8f",conv,k,tol,rate)
                print(file,str)
                # println("Number of projections: $k")
        else
            println("Number of projections: $k")
            @printf("Error = %s\n", tol)
        end
    end
####################################

####################################
    function algorithmMAPmSet(n::Int64,m::Int64,PA,aP,xstar::Vector{Float64},xbar::Vector{Float64},#\file::IOStream,
                                printfile=false,EPS_VAL::Float64=1e-3)
                               # Projecting onto A first
        tol = 1
        k = 0
        conv=String("itmax")
        while (tol >= EPS_VAL)
            # xbar = xstar #Looking for FixingPoint
            for J = 1:m
                xstar = MAP.Projection(PA[J],aP[J],xstar)
            end
            tol = norm(xstar - xbar,2) #True Error

            k+=m
            (k>iter_max) && (warn("Maximum number of projections achieved in MAP"); conv=String("itmax"); break)
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



"""
FindCircumcentermSet(X)

Finds the Circumcenter of vectors ``x_0,x_1,…,x_m``, columns of matrix ``X``,
as described in [^Behling2018a] and [^Behling2018b].

[^Behling2018a]: Behling, R., Bello Cruz, J.Y., Santos, L.-R.: 
Circumcentering the Douglas–Rachford method. Numer. Algorithms. 78(3), 759–776 (2018). 
[doi:10.1007/s11075-017-0399-5](https://doi.org/10.1007/s11075-017-0399-5)
[^Behling2018b]: Behling, R., Bello Cruz, J.Y., Santos, L.-R.: 
On the linear convergence of the circumcentered-reflection method. Oper. Res. Lett. 46(2), 159-162 (2018). 
[doi:10.1016/j.orl.2017.11.018](https://doi.org/10.1016/j.orl.2017.11.018)

"""
    function FindCircumcentermSet(X)
    # Finds the Circumcenter of points X = [X1, X2, X3, ... Xn]
        # println(typeof(X))
        lengthX = length(X)
        if lengthX  == 1
            return X[1]
        elseif lengthX == 2
            return .5*(X[1] + X[2])
        end
        V = []
        b = Float64[]
        # Forms V = [X[2] - X[1] ... X[n]-X[1]]
        # and b = [dot(V[1],V[1]) ... dot(V[n-1],V[n-1])]
        for ind in 2:lengthX
            difXnX1 = X[ind]-X[1]
            push!(V,difXnX1)
            push!(b,dot(difXnX1,difXnX1))
        end

       # Forms Gram Matrix
        dimG = lengthX-1
        G = diagm(b)

        for irow in 1:(dimG-1)
            for icol in  (irow+1):dimG
                G[irow,icol] = dot(V[irow],V[icol])
                G[icol,irow] = G[irow,icol]
            end
        end
        # println(rank(G))
        y = G\b
        # if isposdef(G)
        #     L = cholesky(G)
        #     y = L\b
        # else
        #     @warn "Gram matrix is not SPD"
        #     L = qr(G)
        #     y=L\b
        # end
        CC = X[1]
        for ind in 1:dimG
            CC += .5*y[ind]*V[ind]
        end
        return CC
    end

        ####################################
    """
    ReflectionHyperPlane(x0,a,b)

    Reflects vector ``x_0`` over the hyperplane given by
    ```math
    a⋅x = b.
    ```
    ``a`` must be a column-vector.
    """
    function ReflectionHyperPlane(x0,a,b)
        # Reflection on affine hyperplane a^Tx = b
        alpha = (2.)*((b - dot(x0,a))/dot(a,a))
        return x0 + alpha*a
    end
    ####################################
    """
    BwCRM_Hyperplane(M,b,Blocks; EPS_VAL = 1e-5, xstar = [], sol = [], iter_max = 1e6)

    Applies the CRM method in Blocks as discribed in [^Behling2018a] and [^Behling2019].
    Given `xzero`,  it solves the  *best approximation problem*, i.e., to find the
    closest point `xstar` in ``H`` to `xzero`, where
    ```math
    H = \\bigcap_{i=1}^{m} H_i,
    ```
    and each ``H_i`` is a hyperplane given by `M[i,:]⋅x=b[i].

    In this  version each reflection is computed using the QR factorization.

    [^Behling2018a]	Behling, R., Bello Cruz, J.Y., Santos, L.-R.: On the linear
    convergence of the circumcentered-reflection method. Oper. Res. Lett. 46,
    159–162 (2018). [doi:10.1016/j.orl.2017.11.018](https://doi.org/10.1016/j.orl.2017.11.018)

    [^Behling2019]	Behling, R., Bello Cruz, J.Y., Santos, L.-R.:The Blockwise
    Circumcentered--Reflection Method. (2019) (To be submmited)

    """
    function BwCRM_Hyperplane(M::AbstractArray{Float64,2},b::AbstractArray{Float64,1},
        Blocks::Array{Array{Int64,1},1}; EPS_VAL = 1e-5, xstar = [], sol = [], iter_max = 1e6)
        m,n = size(M)
        NumberBlocks = length(Blocks)
        if isempty(xstar)
            xstar = StartingPoint(n)
        end
        SF = []
        for (index,block) in enumerate(Blocks)
            for row in block
                  F = qr(M[row:row,:]')
                  SF = vcat(SF,F)
             end
        end
        # SF = reshape(SF,:,NumberBlocks)
        normax = 1e2
        iter = 0
        reflec = 0
        while normax > EPS_VAL
            for index in 1:NumberBlocks
                # X = Array{Float64}(undef,n,1)
                X = []
                push!(X,xstar)
                for row in Blocks[index]
                    #Reflection
                    xstar = Reflection(M[row:row,:],b[row:row],SF[row],xstar)
                    push!(X,xstar)
                    reflec +=1
                end
                xstar = FindCircumcentermSet(X)
            end
            iter += 1 #counting outer iterations
            if isempty(sol)
                normax = norm(M*xstar-b,2)
            else
                normax = norm(xstar - sol)
            end
            (iter>iter_max) && (@warn "Maximum number of iterations achieved in BlockWiseCDRM"; conv=String("itmax"); break)
        end
        # @show normax
        return xstar, iter, reflec, normax
    end

        ####################################
    """
    BwCRM_HyperplaneClosedForm(M,b,Blocks; EPS_VAL = 1e-5, xstar = [], sol = [], iter_max = 1e6)

    Applies the CRM method in Blocks as discribed in [^Behling2018a] and [^Behling2019].
    Given `xzero`,  it solves the  *best approximation problem*, i.e., to find
    the closest point `xstar` in ``H`` to `xzero`, where
    ```math
    H = \\bigcap_{i=1}^{m} H_i,
    ```
    and each ``H_i`` is a hyperplane given by `M[i,:]⋅x=b[i].

    In this function,  each reflection is computed using a closed form formula.

    [^Behling2018a]	Behling, R., Bello Cruz, J.Y., Santos, L.-R.: On the linear
    convergence of the circumcentered-reflection method. Oper. Res. Lett. 46,
    159–162 (2018). [doi:10.1016/j.orl.2017.11.018](https://doi.org/10.1016/j.orl.2017.11.018)

    [^Behling2019]	Behling, R., Bello Cruz, J.Y., Santos, L.-R.:The Blockwise
    Circumcentered--Reflection Method. (2019) (To be submmited)

    """
    function BwCRM_HyperplaneClosedForm(M::AbstractArray{Float64,2},b::AbstractArray{Float64,1},
        Blocks::Array{Array{Int64,1},1}; EPS_VAL = 1e-5, xstar = [], sol = [],iter_max = 1e6)
        #Blockwise CRM each subspace is a hyperplane
        m,n = size(M)
        NumberBlocks = length(Blocks)
        if isempty(xstar)
            xstar = StartingPoint(n)
        end
        # xstar = sol + randn(n)
        normax = 1e2
        iter = 0
        reflec = 0
        while normax > EPS_VAL  && iter < iter_max
            for index in 1:NumberBlocks
                X = []
                push!(X,xstar)
                for row in Blocks[index]
                    #Reflection
                    xstar = ReflectionHyperPlane(xstar,M[row:row,:]',b[row])
                    push!(X,xstar)
                    reflec +=1
                end
                xstar = FindCircumcentermSet(X)
            end
            iter += 1 #couting outer iterations
            if isempty(sol)
                normax = norm(M*xstar-b,2)
            else
                normax = norm(xstar - sol)
            end
            (iter>iter_max) && (@warn "Maximum number of iterations achieved in Bw-CDRM"; conv=String("itmax"); break)
        end
        # @show normax
        return xstar, iter, reflec, normax
    end
        ####################################


    """
    MAPmSets(M,b; EPS_VAL = 1e-5, xstar = [], sol = [], iter_max = 1e6)

    Applies the Method of Alternating Projections (MAP) to solve the *best approximation
    problem*, that is, given `xzero`, it finds the closest point `xstar` in ``H`` to `xzero`, where
    ```math
    H = \\bigcap_{i=1}^{m} H_i,
    ```
    and each ``H_i`` is a hyperplane given by `M[i,:]⋅x=b[i].

    In this function,  each reflection is computed using a closed form formula.

    """
    ####################################
    function MAPmSets(M::AbstractArray{Float64,2},b::AbstractArray{Float64,1};
        EPS_VAL = 1e-5, xstar = [], sol = [], iter_max = 1e6)
        #Blockwise CDRM for Hilbert Matrix
        m,n = size(M)
        if isempty(xstar)
            xstar = StartingPoint(n)
        end
        # xstar = sol + randn(n)
        normax = 1e2
        iter = 0
        reflec = 0
        while normax > EPS_VAL
            for row in 1:m
                xstar  += ((b[row] - dot(xstar,M[row,:]))/dot(M[row,:],M[row,:]))*M[row,:] #Reflection
                reflec +=1
            end

            iter += 1 #couting external iterations
            if isempty(sol)
                normax = norm(M*xstar-b,2)
            else
                normax = norm(xstar - sol)
            end
            (iter>iter_max) && (@warn "Maximum number of iterations achieved in BlockWiseCDRM"; conv=String("itmax"); break)
        end
        # @show normax
        return xstar, iter, reflec, normax
    end
    ####################################


    function BwCRM23(M::AbstractArray{Float64,2},b::AbstractArray{Float64,1},
        Blocks::Array{Array{Int64,1},1}; EPS_VAL = 1e-5, xstar=[],sol=[])
        # m,n = size(M)
        # NumberBlocks = length(Blocks)
        if isempty(xstar)
            xstar = StartingPoint(n)
        end
        # xstar = sol + randn(n)
        normax = norm(M*xstar-b,2)
        iter = 0
        reflec = 0
        while normax > EPS_VAL
            X = Array{Float64}(undef,n,1)
            X[:,1] = xstar
            for block in Blocks
                nrowsblock = length(block)
                nsubblocks = div(nrowsblock,2)
                subblocks = formBlocksbyNumber(block[end], nsubblocks, first_index = block[1])
                XB = []
                push!(XB,xstar)
                for rows in subblocks
                    A = Matrix(M[rows,:])
                    SF = qr(A')
                    xstar = Reflection(A,b[rows],SF,xstar)
                    push!(XB,xstar)
                    reflec +=1
                end
                xstar = FindCircumcentermSet(XB)
            end
            iter += 1 #Couting external iterations
            # iter += reflec
            # normax = norm(xstarM-sol,2)
            normax = norm(M*xstar-b,2)
            (iter>iter_max) && (@warn "Maximum number of iterations achieved in BlockWiseCDRM"; conv=String("itmax"); break)
        end
        # @show normax
        return xstar, iter, reflec, normax
    end
    ####################################
    """
    BwCRM(M,b,Blocks; EPS_VAL = 1e-5, xstar = [], sol = [], iter_max = 1e6)

    Applies the CRM method in Blocks as discribed in [^Behling2018a] and [^Behling2019].
    Given `xzero`,  it solves the  *best approximation problem*, i.e., to find
    the closest point `xstar` in ``H`` to `xzero`, where
    ```math
    H = \\bigcap_{i=1}^{p} H_i,
    ```
    and each ``H_i`` is an affine given by rows of matrix M.

    [^Behling2018a]	Behling, R., Bello Cruz, J.Y., Santos, L.-R.: On the linear
    convergence of the circumcentered-reflection method. Oper. Res. Lett. 46,
    159–162 (2018). [doi:10.1016/j.orl.2017.11.018](https://doi.org/10.1016/j.orl.2017.11.018)

    [^Behling2019]	Behling, R., Bello Cruz, J.Y., Santos, L.-R.:The Blockwise
    Circumcentered--Reflection Method. (2019) (To be submmited)
    """
    function BwCRM(M::AbstractArray{Float64,2},b::AbstractArray{Float64,1},
        Blocks::Array{Array{Int64,1},1}; EPS_VAL = 1e-5, xstar=[],sol=[], iter_max = 1e6)
        # If no initial point, compute a random one
        num_row,num_col = size(M)
        if isempty(xstar)
            xstar = StartingPoint(n)
        end
        # Computes the QR factorization of each block
        SF = []
        for block in Blocks
                  F = qr(M[block,:]')
                  SF = vcat(SF,F)
        end
        # xstar = sol + randn(n)
        normax = norm(M*xstar-b,2)
        iter = 0
        reflec = 0
        while normax > EPS_VAL
            X = []
            push!(X,xstar)
            for (index,block) in enumerate(Blocks)
                xstar = Reflection(M[block,:],b[block],SF[index],xstar)
                push!(X,xstar)
                reflec +=1
            end
            xstar = FindCircumcentermSet(X)
            iter += 1 #Couting external iterations
            # iter += reflec
            # normax = norm(xstarM-sol,2)
            normax = norm(M*xstar-b,2)
            (iter>iter_max) && (@warn "Maximum number of iterations achieved in BlockWiseCDRM"; conv=String("itmax"); break)
        end
        # @show normax
        return xstar, iter, reflec, normax
    end
    ####################################

    """
    formBlocksbyNumber(last_index, num_blocks;  first_index = 1)

    Form the indexes of Blocks by number of blocks (`num_blocks`).
    Given the first index `first_index` (default is `1`)  and its last
    index `last_index`, it computes all the Blocks indexes.

    If `num_blocks` is not a divisor of the indexes of the last block
    are complemented until the last index.

    """
    function formBlocksbyNumber(last_index::Int64, num_blocks::Int64;  first_index::Int64 = 1)
        num_index = last_index-first_index + 1
        num_index < num_blocks && error("Number of Blocks exceeds number of elements")
        size_block, add_lastblock = divrem(num_index,num_blocks)
        iszero(add_lastblock) ? status_lastblock = :No : status_lastblock = :Yes
        Blocks = Array{Array{Int64,1}}(undef,0)
        for index in 1:num_blocks
            Blocks = vcat(Blocks,[collect(first_index:first_index+size_block-1)])
            first_index += size_block
        end
        if status_lastblock == :Yes
            Blocks[end] = vcat(Blocks[end],collect(first_index:last_index))
        end
        return Blocks
    end
    ####################################
    """
    formBlocksbySize(last_index, size_block;  first_index = 1)

    Form the indexes of Blocks, by size of each block `size_block`.
    Given the first index `first_index` (default is `1`)  and its last
    index `last_index`, it computes all the Blocks indexes.

    If `size_blocks` is not a divisor of the number of indexes, then the last block
    has the size as the remainder of the division.

    """
    function formBlocksbySize(last_index::Int64, size_block::Int64;  first_index::Int64 = 1)
        num_index = last_index - first_index + 1
        num_blocks, size_lastblock = divrem(num_index,size_block)
        iszero(size_lastblock) ? status_lastblock = :No : status_lastblock = :Yes
        size_block > num_index && error("Size of Block exceeds number of elements")
        Blocks = Array{Array{Int64,1}}(undef,0)
        for index in 1:num_blocks
            Blocks = vcat(Blocks,[collect(first_index:first_index+size_block-1)])
            first_index += size_block
        end
        if status_lastblock == :Yes
            Blocks = vcat(Blocks,[collect(first_index:last_index)])
        end
        return Blocks
    end


    # Makes qr factorization work with adjoint matrices.
    qr(M::Adjoint{Float64,SparseMatrixCSC{Float64,Int64}}) = qr(copy(M))
    qr(M::Adjoint{Float64,Array{Float64,2}}) = qr(copy(M))
# end

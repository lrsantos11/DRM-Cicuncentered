workspace()
include("MAP.jl")
importall MAP

global const EPS_VAL = 1e-3
# fld = @sprintf("../../Experiment_Julia/pairs2/pair%03d.mat",1)
#             include(fld) ## read A,B,cF
#             println(A)
 function GeneralTest(n::Int64 = 100,samples::Int64 = 2,affine::Bool=true)
        # Fix Random
        # srand(2)
        # X = R^n
        for j = 1:samples
            # CC  = GenerateSamples(n,affine)
            # cmax = 3
            # CC  = GenerateRdmPair(n,cmax,affine)
            # cmax = ceil(Integer,n/4)
            fld = @sprintf("../Experiment_Julia/pairs2/pair%03d.mat",j)
            include(fld) ## read A,B,cF
            # println(A)
            ma, n = size(A)
            a = zeros(ma)
            mb, n = size(B)
            b = zeros(mb)
            # A = CC[1]
            # a = CC[2]
            # ma = CC[3]
            # B = CC[4]
            # b= CC[5]
            # mb = CC[6]
            println(rank([A;B]))
            println(size(A))
            println(rank(A))
            println(size(B))
            println(rank(B))
            prob_name  = String("ProblemA$ma""B$mb")
            xzero = randn(n)
            fAngleAB = FriedrichsAngleAB(A,B)
            fAngleAB2 = friedrichs(A,B)
            println(fAngleAB)        
            println(fAngleAB2)
            #########################
            # println("MAP with no acceleration")
            # fname = @sprintf("tables/%d_%.1E-MAP.table",samples,EPS_VAL);
            # open(fname,"a") do file
            #     write(file,prob_name)
            #     time = @elapsed algorithmMAP(A,a,B,b,n,xzero,file,true,1,EPS_VAL)
            #     str = @sprintf(" %10.8e\n",time)
            #     # print(str)
            #     print(file,str)
            # end
            #########################
            # println("MAP with CP acceleration")
            # @time MAP.algorithmMAP(A,a,B,b,n,xzero,2)
            #########################
            # println("MAP with EMAP acceleration")
            # fname = @sprintf("tables/%d_%.1E-EMAP.table",samples,EPS_VAL);
            # open(fname,"a") do file
            #     write(file,prob_name)
            #     time = @elapsed algorithmMAP(A,a,B,b,n,xzero,file,true,3,EPS_VAL)
            #     str = @sprintf(" %10.8e\n",time)
            #     # print(str)
            #     print(file,str)
            # end
            #########################
            # println("Four Point Scheme with no acceleration")
            # fname = @sprintf("tables/%d_%.1E-4pointMAP.table",samples,EPS_VAL);
            # open(fname,"a") do file
            #     write(file,prob_name)
            #     time = @elapsed fourpointscheme(A,a,B,b,n,xzero,file,true,1,EPS_VAL)
            #     str = @sprintf(" %10.8e\n",time)
            #     # print(str)
            #     print(file,str)
            # end
            #########################
            # println("Four Point Scheme with CP acceleration")
            # @time MAP.fourpointscheme(A,a,B,b,n,xzero,2)
            #########################
            # println("Four Point Scheme with EMAP acceleration")
            # fname = @sprintf("tables/%d_%.1E-4pointEMAP.table",samples,EPS_VAL);
            # open(fname,"a") do file
            #     write(file,prob_name)
            #     time = @elapsed fourpointscheme(A,a,B,b,n,xzero,file,true,3,EPS_VAL)
            #     str = @sprintf(" %10.8e\n",time)
            #     # print(str)
            #     print(file,str)
            # end
            #########################
            # println("DRM algorithm")
            # fname = @sprintf("tables/%d_%.1E-DRM.table",samples,EPS_VAL);
            # open(fname,"a") do file
            #     write(file,prob_name)
            #     time = @elapsed algorithmDRM(A,a,B,b,n,xzero,file,true,EPS_VAL)
            #     str = @sprintf(" %10.8e\n",time)
            #     # print(str)
            #     print(file,str)
            # end
            #########################
            # println("DRM Circumcentered algorithm")
            fname = @sprintf("tables/%d_%.1E-DRM-C.table",samples,EPS_VAL);
            open(fname,"a") do file
                write(file,prob_name)
                time = @elapsed algorithmDRM_C(A,a,B,b,n,xzero,file,true,EPS_VAL)
                str = @sprintf(" %10.8e\n",time)
                # print(str)
                print(file,str)
            end
        end
    end

GeneralTest(100,20,false)

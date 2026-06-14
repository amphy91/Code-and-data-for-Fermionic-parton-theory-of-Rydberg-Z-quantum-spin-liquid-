using MKL
using ITensors
using ITensors.HDF5
using Printf
using DelimitedFiles
using LinearAlgebra

function pos(row, col, sl, Lx, Ly)::Int64

    i = 0
    row = (row == 0) ? Ly : ((row == Ly+1) ? 1 : row)
    
    # Snake configuration
    
    if (row < 1 || col < 1 || row > Ly || col > Lx)
        return -1 # indicates that the site is unphysical
    elseif (col%2 == 1)
        i = (col-1)*Ly + row #number of preceding unit cells
        return (i-1)*6+sl
    else
        i = (col-1)*Ly + (Ly-row+1)
        return (i-1)*6+(7-sl)
    end
    
end


function main( )

    Nx = parse(Int64, ARGS[1]) 
    Ny = parse(Int64, ARGS[2]) 
    delta = parse(Int64, ARGS[3]) 
    delta = delta/10.0
    V0 = parse(Int64, ARGS[4]) 
  
    println("$delta $V0")

    N = Nx*Ny*6
    b = floor(Int,Nx/2)*Ny*6
    @show(b)
    
    omega = 1.0 # Rabi frequency
    V_1 = V0 # Blockade interaction
    V_2 = V0 # NNN interaction
    V_3 = V0 # 3NN interaction
    
    sites = siteinds("S=1/2", N; conserve_qns = false)
    
    # Construct the single-site part of the Hamiltonian
    ampo = OpSum()
    for j = 1:N
        ampo += omega,"Sx",j
        ampo += -delta,"projUp",j
    end
    
    # Add a small pinning field
    #double h = 0.001;
    #ampo += -h,"projUp",1; #Favors having a particle on site 1.
  
    # Construct the Rydberg blockade term
    for col = 1:Nx
        for row = 1:Ny
        
            i1 = pos(row,col,1,Nx,Ny)
            i2 = pos(row,col,2,Nx,Ny)
            i3 = pos(row,col,3,Nx,Ny)
            i4 = pos(row,col,4,Nx,Ny)
            i5 = pos(row,col,5,Nx,Ny)
            i6 = pos(row,col,6,Nx,Ny)

            
            ampo += V_1,"projUp",i1, "projUp",i2
            ampo += V_1,"projUp",i2, "projUp",i3
            ampo += V_1,"projUp",i3, "projUp",i1
            ampo += V_1,"projUp",i4, "projUp",i5
            ampo += V_1,"projUp",i5, "projUp",i6
            ampo += V_1,"projUp",i6, "projUp",i4

        end
    end
    
    
    # Construct the second-nearest-neighbor interaction
    for col = 1:2:Nx
        for row = 1:Ny

            i1 = pos(row,col,1,Nx,Ny)
            i2 = pos(row,col,2,Nx,Ny)
            i3 = pos(row,col,3,Nx,Ny)
            i4 = pos(row,col,4,Nx,Ny)
            i5 = pos(row,col,5,Nx,Ny)
            i6 = pos(row,col,6,Nx,Ny)
            
            j1 = pos(row-1,col+1,4,Nx,Ny)
            j3 = pos(row-1,col+1,6,Nx,Ny)
            j5 = pos(row,col+1,1,Nx,Ny)
            j6 = pos(row,col+1,2,Nx,Ny)
            
            if (i2>0 && i4>0)
                ampo += V_2,"projUp",i2, "projUp",i4
            end
            if (i3>0 && i5>0)
                ampo += V_2,"projUp",i3, "projUp",i5
            end
            if (i1>0 && j1>0)
                ampo += V_2,"projUp",i1, "projUp",j1
            end
            if (i3>0 && j3>0)
                ampo += V_2,"projUp",i3, "projUp",j3
            end
            if (i5>0 && j5>0)
                ampo += V_2,"projUp",i5, "projUp",j5
            end
            if (i6>0 && j6>0)
                ampo += V_2,"projUp",i6, "projUp",j6
            end

        end
    end

    for col = 2:2:Nx
        for row = 1:Ny

            i1 = pos(row,col,1,Nx,Ny)
            i2 = pos(row,col,2,Nx,Ny)
            i3 = pos(row,col,3,Nx,Ny)
            i4 = pos(row,col,4,Nx,Ny)
            i5 = pos(row,col,5,Nx,Ny)
            i6 = pos(row,col,6,Nx,Ny)
            
            j1 = pos(row,col+1,4,Nx,Ny)
            j3 = pos(row,col+1,6,Nx,Ny)
            j5 = pos(row+1,col+1,1,Nx,Ny)
            j6 = pos(row+1,col+1,2,Nx,Ny)
            
            if (i2>0 && i4>0)
                ampo += V_2,"projUp",i2, "projUp",i4
            end
            if (i3>0 && i5>0)
                ampo += V_2,"projUp",i3, "projUp",i5
            end
            if (i1>0 && j1>0)
                ampo += V_2,"projUp",i1, "projUp",j1
            end
            if (i3>0 && j3>0)
                ampo += V_2,"projUp",i3, "projUp",j3
            end
            if (i5>0 && j5>0)
                ampo += V_2,"projUp",i5, "projUp",j5
            end
            if (i6>0 && j6>0)
                ampo += V_2,"projUp",i6, "projUp",j6
            end

        end
    end

    
    # Construct the third-nearest-neighbor interaction
    for col = 1:2:Nx
        for row = 1:Ny

            i1 = pos(row,col,1,Nx,Ny)
            i2 = pos(row,col,2,Nx,Ny)
            i3 = pos(row,col,3,Nx,Ny)
            i4 = pos(row,col,4,Nx,Ny)
            i5 = pos(row,col,5,Nx,Ny)
            i6 = pos(row,col,6,Nx,Ny)
            
            j1 = pos(row-1,col+1,6,Nx,Ny)
            j3 = pos(row-1,col+1,4,Nx,Ny)
            j5 = pos(row,col+1,2,Nx,Ny)
            j6 = pos(row,col+1,1,Nx,Ny)
            
            if (i2>0 && i5>0)
                ampo += V_3,"projUp",i2, "projUp",i5
            end
            if (i3>0 && i4>0)
                ampo += V_3,"projUp",i3, "projUp",i4
            end
            if (i1>0 && j1>0)
                ampo += V_3,"projUp",i1, "projUp",j1
            end
            if (i3>0 && j3>0)
                ampo += V_3,"projUp",i3, "projUp",j3
            end
            if (i5>0 && j5>0)
                ampo += V_3,"projUp",i5, "projUp",j5
            end
            if (i6>0 && j6>0)
                ampo += V_3,"projUp",i6, "projUp",j6
            end

        end
    end

    for col = 2:2:Nx
        for row = 1:Ny

            i1 = pos(row,col,1,Nx,Ny)
            i2 = pos(row,col,2,Nx,Ny)
            i3 = pos(row,col,3,Nx,Ny)
            i4 = pos(row,col,4,Nx,Ny)
            i5 = pos(row,col,5,Nx,Ny)
            i6 = pos(row,col,6,Nx,Ny)
            
            j1 = pos(row,col+1,6,Nx,Ny)
            j3 = pos(row,col+1,4,Nx,Ny)
            j5 = pos(row+1,col+1,2,Nx,Ny)
            j6 = pos(row+1,col+1,1,Nx,Ny)
            
            if (i2>0 && i5>0)
                ampo += V_3,"projUp",i2, "projUp",i5
            end
            if (i3>0 && i4>0)
                ampo += V_3,"projUp",i3, "projUp",i4
            end
            if (i1>0 && j1>0)
                ampo += V_3,"projUp",i1, "projUp",j1
            end
            if (i3>0 && j3>0)
                ampo += V_3,"projUp",i3, "projUp",j3
            end
            if (i5>0 && j5>0)
                ampo += V_3,"projUp",i5, "projUp",j5
            end
            if (i6>0 && j6>0)
                ampo += V_3,"projUp",i6, "projUp",j6
            end

        end
    end


    H = MPO(ampo, sites)
    
    # Random state
    psi0 = randomMPS(sites, 10)
    
    
    # Plan to do 200 DMRG sweeps:
    nsweeps = 100;
    # Set maximum MPS bond dimensions for each sweep
    maxdim = [10, 20, 50, 80, 80]
    mindim = [10, 20, 25, 50, 50]
    # Set maximum truncation error allowed when adapting bond dimensions
    cutoff = [1E-12]
    noise = [1E-7, 1E-7, 1E-8, 1E-8, 1E-9, 1E-9, 1E-10, 1E-10, 1E-11, 1E-11, 1E-12, 1E-12, 0.0]
    
    energy, psi = dmrg(H, psi0; nsweeps, maxdim, mindim, cutoff, eigsolve_krylovdim=9)
    psi0 = psi
    #@show energy

    bd = 100
    for nsweeps2 in 1:4
    
        maxdim = bd
        mindim = max(80, floor(Int, bd/2))
        cutoff = [1E-10]

        energy, psi = dmrg(H, psi0; nsweeps=4, maxdim, mindim, cutoff, eigsolve_krylovdim=7)
        psi0 = psi

        orthogonalize!(psi, b)
        U,S,V = svd(psi[b], (linkind(psi, b-1), siteind(psi,b)))
        SvN = 0.0
        for n=1:dim(S, 1)
            p = S[n,n]^2
            SvN -= p * log(p)
        end
        @show(SvN)

        f = h5open(@sprintf("WFZ_%d_%d_D%.1f_V%.1f_B%d.h5", Nx, Ny, delta, V0, bd), "w")
        write(f,"psi",psi)
        close(f)
        
        bd*=2;
        
    end
    
    energy, psi = dmrg(H, psi0; nsweeps=4, maxdim=1600, mindim=800, cutoff=1E-8)
    psi0 = psi

    magz = expect(psi,"Sz")
    open(@sprintf("Sz_%d_%d_D%.1f_V%.1f.txt", Nx, Ny, delta, V0), "w") do io
    for (j,mz) in enumerate(magz)
      println("$j $mz")
      write(io, "$j \t $mz \n")
    end;
    end

    zzcorr = correlation_matrix(psi,"Sz","Sz")
    writedlm(@sprintf("ZZ_%d_%d_D%.1f_V%.1f.txt", Nx, Ny, delta, V0),  zzcorr, ',')

    orthogonalize!(psi, b)
    U,S,V = svd(psi[b], (linkind(psi, b-1), siteind(psi,b)))
    SvN = 0.0
    for n=1:dim(S, 1)
        p = S[n,n]^2
        SvN -= p * log(p)
    end
    @show(SvN)

    open(@sprintf("RubyGS_%d_%d.txt", Nx, Ny), "a") do io
      write(io, "$delta \t $V0 \t $energy \t $SvN \n")
    end;

    # Save wavefunction to disk

    f = h5open(@sprintf("WFZ_%d_%d_D%.1f_V%.1f_B%d.h5", Nx, Ny, delta, V0, bd), "w")
    write(f,"psi",psi)
    close(f)
    
    return nothing
end

main()



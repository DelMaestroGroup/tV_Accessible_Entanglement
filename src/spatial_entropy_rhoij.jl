"""
Calculate both the spatial and the operational entanglement entropies of a
region of consists of sites i and j, using the SVD. The latter is the "entanglement of particles"
introduced by Wiseman and Vaccaro in 2003.
"""

function spatial_entropy_rhoij(basis::AbstractSzbasis, A, d::Vector{Float64}, MaxOccupation::Int)
    B = setdiff(1:basis.K, A)
    # Matrices to SVD
    Amatrices = []

    for i=0:basis.N
        if  i<=MaxOccupation*length(A) && basis.N-i<=MaxOccupation*length(B)
            DimA = num_vectors(basis, i, length(A))
            DimB = num_vectors(basis, basis.N-i, length(B))
            push!(Amatrices, zeros(Float64, DimA, DimB))
        end
    end
    normsDim=min(basis.K-basis.N+1,basis.N+1, MaxOccupation*length(B)+1, MaxOccupation*length(A)+1)
    Index_Shift=0
    if normsDim==min(basis.K-basis.N+1,MaxOccupation*length(B)+1)
       Index_Shift= basis.N -MaxOccupation*length(B)
    end
    norms = zeros(Float64, normsDim )
    for (i, bra) in enumerate(basis)
        AL = zeros(Int, basis.K )
        BL = zeros(Int, basis.K )
        #braA = sub(bra, A)
        #braB = sub(bra, B)
        braA = view(bra, A)
        braB = view(bra, B)
        for k=1: length(A)
            AL[A[k]]=braA[k]
        end
        for k=1: length(B)
            BL[B[k]]=braB[k]
        end
	Flips=0
        SumFlips=0
	for k=1:basis.K 
	    Flips = Flips +(1-2* Flips)* AL[k]
	    SumFlips += Flips* BL[k]
	end
        row = serial_num(basis, length(A), sum(braA), braA)
        col = serial_num(basis, length(B), sum(braB), braB)
        Amatrices[1 + sum(braA)-Index_Shift][row, col] = d[i]*(-1)^ SumFlips
        norms[1 + sum(braA)-Index_Shift] += d[i]^2
    end

    norm_err = abs(sum(norms) - 1.0)
    if norm_err > 1e-12
        warn("norm error ", norm_err)
    end

    rho11=Float64
    rho44=Float64
    rho22=Float64
    rho33=Float64
    rho23=Float64
    rho32=Float64

    for (i, Amatrix) in enumerate(Amatrices)
      rho_n= Amatrix *transpose(Amatrix)
      if i==1
        rho11= rho_n[1,1]
      elseif i==3
        rho44= rho_n[1,1]
      else
        rho22= rho_n[2,2]
        rho33= rho_n[1,1]
        rho23= rho_n[2,1]
        rho32= rho_n[1,2]
      end 
    end
    rho23_err = abs(rho23-rho32)
    if rho23_err > 1e-12
        warn("rho23 error ", rho23_err)
    end
    rho11,rho22,rho33,rho44,rho23
end
spatial_entropy_rhoij(basis::AbstractSzbasis, sitei::Int, sitej::Int, d::Vector{Float64}, MaxOccupation::Int) = spatial_entropy_rhoij(basis, [sitei,sitej], d, MaxOccupation)

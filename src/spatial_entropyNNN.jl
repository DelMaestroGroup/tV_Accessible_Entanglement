"""
Calculate both the spatial and the operational entanglement entropies of a
region A, using the SVD. The latter is the "entanglement of particles"
introduced by Wiseman and Vaccaro in 2003.
"""
function spatial_entropyNNN(basis::AbstractSzbasis, A, d::Vector{Float64}, MaxOccupation::Int)
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
        Amatrices[1 + sum(braA)-Index_Shift][row, col] = d[i]
        norms[1 + sum(braA)-Index_Shift] += d[i]^2
    end

    norm_err = abs(sum(norms) - 1.0)

    if norm_err > 1e-12
        warn("norm error ", norm_err)
    end
    Ss_raw = [svdvals!(Amatrix) for Amatrix in Amatrices]


    # Spatial.
    S_sp = vcat(Ss_raw...)

    err_sp = abs(sum(S_sp.^2) - 1.0)

    if err_sp > 1e-12
        warn("RDM eigenvalue error ", err_sp)
    end

    S2_sp = -log(sum(S_sp.^4))
    S1_sp=0
    for k=1:length(S_sp)
        if S_sp[k]^2>0
            S1_sp -=S_sp[k]^2*log(S_sp[k]^2)
	end
    end
    # Operational.
    Ss_op = [S / sqrt(n) for (S, n) in zip(Ss_raw, norms)]

    errs_op = [abs(sum(S.^2) - 1.0) for S in Ss_op]

    if any(errs_op .> 1e-12)
        warn("RDM eigenvalue error_OP ", maximum(errs_op))
    end

    S2s_op = [-log(sum(S.^4)) for S in Ss_op]
    S2_op = dot(norms, S2s_op)
    S1s_op = [0.0 for S in Ss_op]
    for (i, S) in enumerate(Ss_op)
       for k=1:length(S)
          if S[k]^2>0
              S1s_op[i] -=S[k]^2*log(S[k]^2)
	  end
       end
    end
    S1_op = dot(norms, S1s_op)
    N2=0
    N1=0
    for i=1:normsDim
       N1+=(i-1+ Index_Shift)*norms[i]
       N2+=(i-1+ Index_Shift)^2*norms[i]
    end
    Sigma2_n=N2-N1^2
    S1_sp, S1_op, S2_sp, S2_op,Sigma2_n


end
spatial_entropyNNN(basis::AbstractSzbasis, Asize::Int, d::Vector{Float64}, MaxOccupation::Int) = spatial_entropyNNN(basis, 1:2:2*Asize, d, MaxOccupation)

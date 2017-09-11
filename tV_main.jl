#Renyi entanglement entropy of Bose-Hubbard chains in 1D.

push!(LOAD_PATH, joinpath(dirname(@__FILE__), "src"))

using tVDiagonalize
using ArgParse
using JeszenszkiBasis

s = ArgParseSettings()
s.autofix_names = true
@add_arg_table s begin
    "M"
        help = "number of sites"
        arg_type = Int
        required = true
    "N"
        help = "number of particles"
        arg_type = Int
        required = true
    "--out"
        metavar = "FILE"
        help = "path to output file"
        required = true
    "--site-max"
        metavar = "N"
        help = "site occupation restriction"
        arg_type = Int
end

add_arg_group(s, "boundary conditions")
@add_arg_table s begin
    "--pbc"
        help = "periodic boundary conditions (default)"
        arg_type = BdryCond
        action = :store_const
        dest_name = "boundary"
        constant = PBC
        default = PBC
    "--obc"
        help = "open boundary conditions"
        arg_type = BdryCond
        action = :store_const
        dest_name = "boundary"
        constant = OBC
        default = PBC
end

add_arg_group(s, "two sites")
@add_arg_table s begin
    "--notsrdm"
        help = "do not generate the two sites RDM"
        arg_type = twositesRDM
        action = :store_const
        dest_name = "twosites"
        constant = NOTSRDM
        default = NOTSRDM
    "--tsrdm"
        help = "generate the two sites RDM"
        arg_type = twositesRDM
        action = :store_const
        dest_name = "twosites"
        constant = TSRDM
        default = NOTSRDM
end

add_arg_group(s, "Next nearest neighbor bipartition")
@add_arg_table s begin
    "--nonnn"
        help = "do not include next nearest neighbor bipartition"
        arg_type = NNNbipart
        action = :store_const
        dest_name = "nnnbip"
        constant = NONNN
        default = NONNN
    "--nnn"
        help = "include next nearest neighbor bipartition"
        arg_type = NNNbipart
        action = :store_const
        dest_name = "nnnbip"
        constant = NNN
        default = NONNN
end



add_arg_group(s, "BH parameters")
@add_arg_table s begin
    "--u-min"
        metavar = "U"
        help = "minimum U"
        arg_type = Float64
        default = 1.0
    "--u-max"
        metavar = "U"
        help = "maximum U"
        arg_type = Float64
        default = 20.0
    "--u-step"
        metavar = "U"
        help = "U step"
        arg_type = Float64
    "--u-num"
        metavar = "N"
        help = "number of U"
        arg_type = Int
    "--u-log"
        help = "use logarithmic scale for U"
        action = :store_true
    "--t"
        metavar = "t"
        help = "t value"
        arg_type = Float64
        default = 1.0
    "--j-max"
        metavar = "j_max"
        help = "maximum j"
        arg_type = Int
        default = -1
end 

add_arg_group(s, "entanglement entropy")
@add_arg_table s begin
    "--ee"
        metavar = "XA"
        help = "compute all EEs with partition size XA"
        arg_type = Int
        required = true
end
c = parsed_args = parse_args(ARGS, s, as_symbols=true)

# Number of sites
const M = c[:M]
# Number of particles
const N = c[:N]
# Output file
const output = c[:out]
# Site occupation restriction
const site_max = c[:site_max]
# Boundary conditions
const boundary = c[:boundary]
# Size of region A
const Asize = c[:ee]

const twosites = c[:twosites]
const nnnbip = c[:nnnbip]
if nnnbip==NNN && Asize>M/2
    println("for --nnn, --ee must be <= M/2")
    exit(1)
end
if c[:j_max]==-1
   r_max=M-1
else
   r_max=c[:j_max]-1
end
if c[:u_log] && c[:u_num] === nothing
    println("--u-log must be used with --u-num")
    exit(1)
end

if c[:u_step] === nothing
    if c[:u_num] === nothing
        U_range = c[:u_min]:0.5:c[:u_max]
    else
        if c[:u_log]
            U_range = logspace(c[:u_min], c[:u_max], c[:u_num])
        else
            U_range = linspace(c[:u_min], c[:u_max], c[:u_num])
        end
    end
else
    if c[:u_num] === nothing
        U_range = c[:u_min]:c[:u_step]:c[:u_max]
    else
        println("--u-step and --u-num may not both be supplied")
        exit(1)
    end
end

if site_max === nothing
    const basis = Szbasis(M, N)
else
    const basis = RestrictedSzbasis(M, N, site_max)
end

#_______________________
Afiles = []

ll=length(basis)
rho11 =zeros(Float64, r_max)
rho22 =zeros(Float64, r_max)
rho33 =zeros(Float64, r_max)
rho44 =zeros(Float64, r_max)
rho23 =zeros(Float64, r_max)

Vz=zeros(Float64, ll)
Vh=zeros(Float64, ll)
Vl=zeros(Float64, ll)
wf=zeros(Float64, ll)

for i=1:ll
   Vh[i]=0.1
   #Vl[i]=0.1
   Vl[i]=0.0
   Vz[i]=1.0
end

if boundary==OBC
    num_links= basis.K-1
elseif boundary==PBC
    num_links= basis.K
end

for (i, bra) in enumerate(basis)
    cc=0
    for j=1: num_links
        j_next = j % basis.K + 1
        cc+=bra[j]*bra[j_next]
    end
    if cc== basis.N-1
        Vl[serial_num(basis, bra)]=1.0
    elseif cc==0
        Vh[serial_num(basis, bra)]=1.0
    end
end

#Normh=0
Norml=0
#Normz=0
for i=1:ll
#   Normh+=Vh[i]^2
   Norml+=Vl[i]^2
#   Normz+=Vz[i]^2
end

for i=1:ll
# Vh[i]/=sqrt(Normh)
 Vl[i]/=sqrt(Norml)
# Vz[i]/=sqrt(Normz)
end

#_______________________
if twosites== TSRDM
   output11= output[1:end-4]"_r11"output[end-3:end]
   f11=open(output11, "w") 
        write(f11, "# M=$(M), N=$(N), max=$(site_max), $(boundary), r_max=$(r_max)\n")
        write(f11, "# U/t E0/t    ")
        for r=1:r_max
            write(f11, " rho11($(r))    ")
        end  
        write(f11, "\n")

   output22= output[1:end-4]"_r22"output[end-3:end]
   f22=open(output22, "w") 
        write(f22, "# M=$(M), N=$(N), max=$(site_max), $(boundary), r_max=$(r_max)\n")
        write(f22, "# U/t E0/t    ")
        for r=1:r_max
            write(f22, " rho22($(r))    ")
        end  
        write(f22, "\n")

   output33= output[1:end-4]"_r33"output[end-3:end]
   f33=open(output33, "w") 
        write(f33, "# M=$(M), N=$(N), max=$(site_max), $(boundary), r_max=$(r_max)\n")
        write(f33, "# U/t E0/t    ")
        for r=1:r_max
            write(f33, " rho33($(r))    ")
        end  
        write(f33, "\n")

   output44= output[1:end-4]"_r44"output[end-3:end]
   f44=open(output44, "w") 
        write(f44, "# M=$(M), N=$(N), max=$(site_max), $(boundary), r_max=$(r_max)\n")
        write(f44, "# U/t E0/t    ")
        for r=1:r_max
            write(f44, " rho44($(r))    ")
        end  
        write(f44, "\n")

   output23= output[1:end-4]"_r23"output[end-3:end]
   f23=open(output23, "w") 
        write(f23, "# M=$(M), N=$(N), max=$(site_max), $(boundary), r_max=$(r_max)\n")
        write(f23, "# U/t E0/t    ")
        for r=1:r_max
            write(f23, " rho23($(r))    ")
        end  
        write(f23, "\n")
end
if nnnbip== NNN
   outputNNN= output[1:end-4]"NNN"output[end-3:end]
   fNNN=open(outputNNN, "w") 
   if site_max === nothing
       write(fNNN, "# M=$(M), N=$(N), $(boundary)\n")
   else
       write(fNNN, "# M=$(M), N=$(N), max=$(site_max), $(boundary),   Next nearest neighbor\n")
   end
   write(fNNN, "# U/t E0/t              S1(l=$(Asize))            S1_OP(l=$(Asize))          S2(l=$(Asize))            S2_OP(l=$(Asize))              Sigma2_n(l=$(Asize))\n")
end




open(output, "w") do f
    if site_max === nothing
        write(f, "# M=$(M), N=$(N), $(boundary)\n")
    else
        write(f, "# M=$(M), N=$(N), max=$(site_max), $(boundary)\n")


    end
    write(f, "# U/t E0/t              S1(l=$(Asize))            S1_OP(l=$(Asize))          S2(l=$(Asize))            S2_OP(l=$(Asize))              Sigma2_n(l=$(Asize))\n")


           # wf=Vh
            wf=Vl
           # wf=Vz

    for U in U_range

        # Create the Hamiltonian
        H = sparse_hamiltonian(basis, c[:t], U, boundary=boundary)
#	print(" sparse_hamiltonian finish\n ")

        # Perform the Lanczos diagonalization to obtain the lowest eigenvector
        # http://docs.julialang.org/en/release-0.3/stdlib/linalg/?highlight=lanczos
        d = eigs(H, nev=1, which=:SR,tol=1e-12,v0=wf)
        wf = vec(d[2][1:ll])

        # Calculate the second Renyi entropy
        s1_spatial, s1_operational, s2_spatial, s2_operational, Sigma2_n = spatial_entropy(basis, Asize, wf, site_max)

        write(f, "$(U/c[:t]) $(d[1][1]/c[:t]) $(s1_spatial) $(s1_operational) $(s2_spatial) $(s2_operational) $(Sigma2_n)\n")
        flush(f)
        if nnnbip== NNN
           s1_spatial, s1_operational, s2_spatial, s2_operational, Sigma2_n = spatial_entropyNNN(basis, Asize, wf, site_max)
           write(fNNN, "$(U/c[:t]) $(d[1][1]/c[:t]) $(s1_spatial) $(s1_operational) $(s2_spatial) $(s2_operational) $(Sigma2_n)\n")
           flush(fNNN)
        end

        if twosites== TSRDM
           for r=1:r_max
              rho11[r],rho22[r],rho33[r],rho44[r],rho23[r] = spatial_entropy_rhoij(basis, 1,r+1, wf, site_max)
           end

           write(f11, "$(U/c[:t]) $(d[1][1]/c[:t])   ")              
           for r=1:r_max
               write(f11, " $(rho11[r])    ")
           end  
           write(f11, "\n")
           flush(f11)

           write(f22, "$(U/c[:t]) $(d[1][1]/c[:t])   ")              
           for r=1:r_max
               write(f22, " $(rho22[r])    ")
           end  
           write(f22, "\n")
           flush(f22)

           write(f33, "$(U/c[:t]) $(d[1][1]/c[:t])   ")              
           for r=1:r_max
               write(f33, " $(rho33[r])    ")
           end  
           write(f33, "\n")
           flush(f33)

           write(f44, "$(U/c[:t]) $(d[1][1]/c[:t])   ")              
           for r=1:r_max
               write(f44, " $(rho44[r])    ")
           end  
           write(f44, "\n")
           flush(f44)

           write(f23, "$(U/c[:t]) $(d[1][1]/c[:t])   ")              
           for r=1:r_max
               write(f23, " $(rho23[r])    ")
           end  
           write(f23, "\n")
           flush(f23)

        end
    end 
   if twosites== TSRDM
      close(f11)
      close(f22)
      close(f33)
      close(f44)
      close(f23)
   end
   if nnnbip== NNN
      close(fNNN)
   end
 end

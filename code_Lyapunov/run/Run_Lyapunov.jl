#=
julia --check-bounds=no run/Run_Lyapunov.jl \
--Mtot 1.0 --Nv 20 --Npart 100 \
--kmin 10 --kmax 100 \
--vmax 6.0 --dt 1 \
--Nsteps 1000 --hierarchy F \
--IC MAXW --method FFT
=#

##################################################
include("../src/Main.jl") # Loading the source code
using Plots               # For plotting
##################################################


##################################################
# Initial dispersions (anisotropic Maxwellian)
##################################################
const T_X = 1.5 * T_0
const T_Y = 0.5 * T_0

##################################################
# Write the results on a hf5 file
##################################################

function write_file_Lyapunov_eq()
    #####
    namefile = (@__DIR__) * "/../data/" *
                "data_Lyapunov_HIERARCHY_" * string(HIERARCHY) *
                    "_NPART_" * string(NPART) * 
                    "_NV_" * string(NV) * 
                    "_IC_" * string(IC) * 
                    "_METHOD_" * string(METHOD) * 
               ".hf5"

    #####
    file = h5open(namefile,"w") # Opening the file
    #####
    # Writing the parameters of the run
    write(file,"HIERARCHY", HIERARCHY)
    write(file,"METHOD", METHOD)
    write(file,"IC", IC)
    write(file,"NV", NV)
    write(file,"DV", DV)
    write(file,"VMAX", VMAX)
    write(file,"DT", DT)
    write(file,"NSTEPS", NSTEPS)  
    write(file,"KMAX", KMAX)
    write(file,"KMIN", KMIN)
    write(file,"G", G)
    write(file,"NPART", NPART)
    write(file,"MASS_PART", MASS_PART)
    #####
    # Writing the data of the run
    write(file,"TAB_TIME", TAB_TIME)
    write(file,"TAB_TXX_AVG", TAB_Tx)
    write(file,"TAB_TYY_AVG", TAB_Ty)
    #
    # write(file,"TAB_OBSERV", TAB_OBSERV)
    #
    write(file,"TAB_SXXX_AVG", TAB_SXXX)
    write(file,"TAB_SXXY_AVG", TAB_SXXY)
    write(file,"TAB_SXYY_AVG", TAB_SXYY)
    write(file,"TAB_SYYY_AVG", TAB_SYYY)
    #
    write(file,"TAB_NONPOLY_AVG", TAB_NON_POLY_OBS)
    #
    write(file,"TAB_QXXXX_AVG", TAB_QXXXX)
    write(file,"TAB_QXXXY_AVG", TAB_QXXXY)
    write(file,"TAB_QXXYY_AVG", TAB_QXXYY)
    write(file,"TAB_QXYYY_AVG", TAB_QXYYY)
    write(file,"TAB_QYYYY_AVG", TAB_QYYYY)
    #
    write(file,"TAB_POLYSINCOS_AVG", TAB_POLYSINCOS)
    write(file,"TAB_GAUSSIAN_AVG", TAB_GAUSSIAN)
    #
    write(file,"TAB_COSPLUS_AVG", TAB_COS_PLUS)
    write(file,"TAB_COSMINUS_AVG", TAB_COS_MINUS)
    write(file,"TAB_VSINVSIN_AVG", TAB_VSIN_VSIN)
    write(file,"TAB_TANHVX_AVG", TAB_TANH_VX)
    write(file,"TAB_SIGMOIDVX_AVG", TAB_SIGMOID_VX)
    write(file,"TAB_TANHANISO_AVG", TAB_TANH_ANISO)
    #
    if HIERARCHY == "S2"
        write(file,"TAB_VAR_M", TAB_VAR_M)
        write(file,"TAB_VAR_E", TAB_VAR_E)
        write(file,"TAB_VAR_Px", TAB_VAR_Px)
        write(file,"TAB_VAR_Py", TAB_VAR_Py)
        write(file,"TAB_VAR_Ptot", TAB_VAR_Ptot)
        #
        write(file,"TAB_TXX_VAR", TAB_VAR_Tx)
        write(file,"TAB_TYY_VAR", TAB_VAR_Ty)
        #
        write(file,"TAB_SXXX_VAR", TAB_VAR_SXXX)
        write(file,"TAB_SXXY_VAR", TAB_VAR_SXXY)
        write(file,"TAB_SXYY_VAR", TAB_VAR_SXYY)
        write(file,"TAB_SYYY_VAR", TAB_VAR_SYYY)
        #
        write(file,"TAB_NONPOLY_VAR", TAB_VAR_NON_POLY_OBS)
        #
        write(file,"TAB_QXXXX_VAR", TAB_VAR_QXXXX)
        write(file,"TAB_QXXXY_VAR", TAB_VAR_QXXXY)
        write(file,"TAB_QXXYY_VAR", TAB_VAR_QXXYY)
        write(file,"TAB_QXYYY_VAR", TAB_VAR_QXYYY)
        write(file,"TAB_QYYYY_VAR", TAB_VAR_QYYYY)
        #
        write(file,"TAB_POLYSINCOS_VAR", TAB_VAR_POLYSINCOS)
        write(file,"TAB_GAUSSIAN_VAR", TAB_VAR_GAUSSIAN)
        #
        write(file,"TAB_COSPLUS_VAR", TAB_VAR_COS_PLUS)
        write(file,"TAB_COSMINUS_VAR", TAB_VAR_COS_MINUS)
        write(file,"TAB_VSINVSIN_VAR", TAB_VAR_VSIN_VSIN)
        write(file,"TAB_TANHVX_VAR", TAB_VAR_TANH_VX)
        write(file,"TAB_SIGMOIDVX_VAR", TAB_VAR_SIGMOID_VX)
        write(file,"TAB_TANHANISO_VAR", TAB_VAR_TANH_ANISO)
    end
    #
    if HIERARCHY == "S3"
        write(file,"TAB_M3_M", TAB_M3_M)
        write(file,"TAB_M3_E", TAB_M3_E)
        write(file,"TAB_M3_Px", TAB_M3_Px)
        write(file,"TAB_M3_Py", TAB_M3_Py)
        write(file,"TAB_M3_Ptot", TAB_M3_Ptot)
        #
        write(file,"TAB_TXX_M3", TAB_M3_Tx)
        write(file,"TAB_TYY_M3", TAB_M3_Ty)
        #
        write(file,"TAB_SXXX_M3", TAB_M3_SXXX)
        write(file,"TAB_SXXY_M3", TAB_M3_SXXY)
        write(file,"TAB_SXYY_M3", TAB_M3_SXYY)
        write(file,"TAB_SYYY_M3", TAB_M3_SYYY)
        #
        write(file,"TAB_NONPOLY_M3", TAB_M3_NON_POLY_OBS)
        #
        write(file,"TAB_QXXXX_M3", TAB_M3_QXXXX)
        write(file,"TAB_QXXXY_M3", TAB_M3_QXXXY)
        write(file,"TAB_QXXYY_M3", TAB_M3_QXXYY)
        write(file,"TAB_QXYYY_M3", TAB_M3_QXYYY)
        write(file,"TAB_QYYYY_M3", TAB_M3_QYYYY)
        #
        write(file,"TAB_POLYSINCOS_M3", TAB_M3_POLYSINCOS)
        write(file,"TAB_GAUSSIAN_M3", TAB_M3_GAUSSIAN)
        #
        write(file,"TAB_COSPLUS_M3", TAB_M3_COS_PLUS)
        write(file,"TAB_COSMINUS_M3", TAB_M3_COS_MINUS)
        write(file,"TAB_VSINVSIN_M3", TAB_M3_VSIN_VSIN)
        write(file,"TAB_TANHVX_M3", TAB_M3_TANH_VX)
        write(file,"TAB_SIGMOIDVX_M3", TAB_M3_SIGMOID_VX)
        write(file,"TAB_TANHANISO_M3", TAB_M3_TANH_ANISO)
    end
    #####
    # Writing the invariants
    write(file,"TAB_E", TAB_E)
    write(file,"TAB_M", TAB_M)
    write(file,"TAB_Px", TAB_Px)
    write(file,"TAB_Py", TAB_Py)
    #####  
    close(file) # Output
end


##################################################
# Plot of Txx: mean, variance (S2) or M3 (S3)
##################################################

function plot_Txx()
    #####
    # Choice of the quantity to plot depending on HIERARCHY
    if HIERARCHY == "S2"
        tab    = TAB_VAR_Tx
        label  = "Var(Tx)"
        ylabel = "Variance of Tx"
        title  = "Variance of Tx"
        tag    = "VAR"
    elseif HIERARCHY == "S3"
        tab    = TAB_M3_Tx
        label  = "M3(Tx)"
        ylabel = "Third moment of Tx"
        title  = "Third moment M3 of Tx"
        tag    = "M3"
    else # Hierarchy F: mean
        tab    = TAB_Tx
        label  = "Tx"
        ylabel = "Temperature"
        title  = "Temperature relaxation Tx"
        tag    = "AVG"
    end
    #####
    plt = plot(TAB_TIME, tab,
               label=label,
               xlabel="Time (t)",
               ylabel=ylabel,
               title=title,
               linewidth=2,
               grid=true)
    #####
    namefig = (@__DIR__) * "/../data/" *
              "plot_TXX_" * tag *
              "_HIERARCHY_" * string(HIERARCHY) *
              "_NPART_" * string(NPART) *
              "_NV_" * string(NV) *
              "_IC_" * string(IC) *
              "_METHOD_" * string(METHOD) *
              ".png"
    savefig(plt, namefig)
    println("Figure saved | ", namefig)
    #####
    return plt
end


##################################################
# Performing the run
##################################################
function run_simulation!()
    # Initialization
    init!(T_X, T_Y)
    TAB_TIME[1] = 0.0
    store_observables!(1)
    println("Initial mass = ", TAB_M[1])

    # Time stepping loop
    for n in 1:NSTEPS
        integrate_dt!() # integrate S over [t, t+dt]
        TAB_TIME[n+1] = n * DT
        println(n)
        store_observables!(n+1)
    end
end

##################################################
# Running the code
##################################################
print("Threads used |")
println(Threads.nthreads())
print("Integrating |")
@time run_simulation!()
write_file_Lyapunov_eq()
plot_Txx()
#
println("Final time | ",TAB_TIME[NSTEPS+1]," | err_Ptot | ",compute_Ptot()," | err_Etot | ",compute_err_Etot())
if HIERARCHY == "S3"
    println("err_M3_mass | ", compute_err_M3_mass()," | err_M3_Ekin | ", compute_err_M3_Etot()   )
    println("MEMORY TOT S3 | ",Base.format_bytes(4*Base.summarysize(TAB_S3)))
elseif HIERARCHY == "S2"
    println("err_M2_mass | ", compute_err_var_mass()," | err_M2_Ekin | ", compute_err_var_Etot()   )
    println("MEMORY TOT S2 | ",Base.format_bytes(10*Base.summarysize(TAB_S2)))
end

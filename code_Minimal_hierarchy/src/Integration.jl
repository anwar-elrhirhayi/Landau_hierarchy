##################################################
# Integrates for one timestep
##################################################

if HIERARCHY == "F"
    ####
    if METHOD == "FD"
        include("../src/src_Hierarchy_F/Landau_FD.jl")
    elseif METHOD == "FFT"
        include("../src/src_Hierarchy_F/Landau_FFT.jl")
    end
    ####
    include("../src/src_Hierarchy_F/Initialisation_F.jl")
    include("../src/src_Hierarchy_F/Integration_F.jl")

    const init! = init_Hierarchy_F!
    const integrate_dt! = integrate_dt_hierarchy_F!
    ####
elseif HIERARCHY == "S2"
    ####
    include("../src/src_Hierarchy_S2/Arrays_S2.jl")
    include("../src/src_Hierarchy_S2/Initialisation_S2.jl")
    include("../src/src_Hierarchy_S2/Gradients_S2.jl")
    include("../src/src_Hierarchy_S2/Integration_S2.jl")
    include("../src/src_Hierarchy_S2/Landau_correction_S2.jl")
    if METHOD == "FD"
        include("../src/src_Hierarchy_S2/Lyapunov_S2_FD.jl")
    elseif METHOD == "FFT"
        include("../src/src_Hierarchy_S2/Lyapunov_S2_FFT.jl")
    end 
    include("../src/src_Hierarchy_S2/Source_S2.jl")

    const init! = init_Hierarchy_S2!
    const integrate_dt! = integrate_dt_Hierarchy_S2!
    ####
elseif HIERARCHY == "S3"
    ####
    include("../src/src_Hierarchy_S3/Arrays_S3.jl")
    include("../src/src_Hierarchy_S3/Initialisation_S3.jl")
    include("../src/src_Hierarchy_S3/Landau_S3.jl")
    include("../src/src_Hierarchy_S3/Lyapunov_S2.jl")
    ####
    if METHOD == "FD"
        include("../src/src_Hierarchy_S3/Lyapunov_S3_FD.jl")
        include("../src/src_Hierarchy_S3/Integrals_S3_FD.jl")    
    elseif METHOD == "FFT"
        include("../src/src_Hierarchy_S3/Lyapunov_S3_FFT.jl")
        include("../src/src_Hierarchy_S3/Integrals_S3_FFT.jl")    
    end 
    ####   
    include("../src/src_Hierarchy_S3/Gradients_S2.jl")
    include("../src/src_Hierarchy_S3/NL_Lyapunov_S3.jl")
    include("../src/src_Hierarchy_S3/Source_S2.jl")
    include("../src/src_Hierarchy_S3/Source_S3.jl")
    include("../src/src_Hierarchy_S3/Integration_S3_symmetry.jl")

    const init! = init_Hierarchy_S3!
    const integrate_dt! = integrate_dt_Hierarchy_S3!
    ####
end  
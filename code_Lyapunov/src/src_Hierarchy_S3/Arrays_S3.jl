##################################################
# PREALLOCATED ARRAYS FOR PRECOMPUTE INTEGRALS
##################################################

const TAB_INT_Bxx_F   = zeros(Float64, NV2)
const TAB_INT_Bxy_F   = zeros(Float64, NV2)
const TAB_INT_Byy_F   = zeros(Float64, NV2)
#
const TAB_INT_Bxx_dFx = zeros(Float64, NV2)
const TAB_INT_Bxy_dFx = zeros(Float64, NV2)
const TAB_INT_Bxy_dFy = zeros(Float64, NV2)
const TAB_INT_Byy_dFy = zeros(Float64, NV2)


##################################################
# PREALLOCATED ARRAYS FOR S2
##################################################

##################################################
# Lyapunov Linear S2 part
##################################################

# ATTENTION, here TAB_S2 correspond to S_2^(1) and not 1/N S_2^(1)
const TAB_S2            = zeros(Float64, NV2, NV2)    # 2 point correlation fonction
const TAB_S2_SLICE      = zeros(Float64, NV2)         # S(v, v' fixed) 
const TAB_dS2_SLICE_dx  = zeros(Float64, NV2)               
const TAB_dS2_SLICE_dy  = zeros(Float64, NV2)         
const TAB_dS2_dt        = zeros(Float64, NV2, NV2)   

const TAB_dS2_dx        = zeros(Float64, NV2, NV2)              
const TAB_dS2_dy        = zeros(Float64, NV2, NV2)   


##################################################
# Lyapunov Source S2 part 
##################################################

const TAB_SOURCE_S2 = zeros(Float64, NV2, NV2) 

##################################################
# Observables and Invariants
##################################################
const TAB_M3_M   = zeros(Float64, NSTEPS+1) # M3 of total system mass
const TAB_M3_E   = zeros(Float64, NSTEPS+1) # M3 of total system kinetic energy
const TAB_M3_Px  = zeros(Float64, NSTEPS+1)   # M3 of total momentum x-direction
const TAB_M3_Py  = zeros(Float64, NSTEPS+1)   # M3 of total momentum y-direction
const TAB_M3_Ptot  = zeros(Float64, NSTEPS+1) # M3 of total momentum

# Temperature
const TAB_M3_Tx  = zeros(Float64, NSTEPS+1)  # M3 of Tx
const TAB_M3_Ty  = zeros(Float64, NSTEPS+1)  # M3 of Ty

# Cubic Observable
const TAB_M3_SXXX = zeros(Float64,NSTEPS+1)
const TAB_M3_SXXY = zeros(Float64,NSTEPS+1)
const TAB_M3_SXYY = zeros(Float64,NSTEPS+1)
const TAB_M3_SYYY = zeros(Float64,NSTEPS+1)
# Non polynomial observables
const TAB_M3_NON_POLY_OBS = zeros(Float64,NSTEPS+1)

# Quartic Observable M3 arrays
const TAB_M3_QXXXX = zeros(Float64,NSTEPS+1)
const TAB_M3_QXXXY = zeros(Float64,NSTEPS+1)
const TAB_M3_QXXYY = zeros(Float64,NSTEPS+1)
const TAB_M3_QXYYY = zeros(Float64,NSTEPS+1)
const TAB_M3_QYYYY = zeros(Float64,NSTEPS+1)
# Poly-sin-cos and Gaussian observables M3 arrays
const TAB_M3_POLYSINCOS = zeros(Float64,NSTEPS+1)
const TAB_M3_GAUSSIAN   = zeros(Float64,NSTEPS+1)
# Anisotropic cos channel pair, vsin_vsin, tanh/sigmoid observables M3 arrays
const TAB_M3_COS_PLUS    = zeros(Float64,NSTEPS+1)
const TAB_M3_COS_MINUS   = zeros(Float64,NSTEPS+1)
const TAB_M3_VSIN_VSIN   = zeros(Float64,NSTEPS+1)
const TAB_M3_TANH_VX     = zeros(Float64,NSTEPS+1)
const TAB_M3_SIGMOID_VX  = zeros(Float64,NSTEPS+1)
const TAB_M3_TANH_ANISO  = zeros(Float64,NSTEPS+1)

##################################################
# PREALLOCATED ARRAYS for S3
##################################################

# ATTENTION, here TAB_S3 correspond to S_3^(2) and not 1/N S_3^(2)
const TAB_S3            = zeros(Float64, NV2, NV2, NV2)    # 3 point correlation fonction 
const TAB_dS3_dt        = zeros(Float64, NV2, NV2, NV2)   

const TAB_S3_SLICE      = zeros(Float64, NV2)         
const TAB_dS3_SLICE_dx  = zeros(Float64, NV2)               
const TAB_dS3_SLICE_dy  = zeros(Float64, NV2)  

const TAB_S3_SLICE_PARALLEL = [zeros(Float64, NV2) for tid = 1:Threads.nthreads()]
const TAB_DS3_DX_PARALLEL   = [zeros(Float64, NV2) for tid = 1:Threads.nthreads()]
const TAB_DS3_DY_PARALLEL   = [zeros(Float64, NV2) for tid = 1:Threads.nthreads()]
##################################################
# Lyapunov S3 Linear part
##################################################       

const TAB_dS3_dx1       = zeros(Float64, NV2, NV2)              
const TAB_dS3_dy1       = zeros(Float64, NV2, NV2)   
const TAB_dS3_dx2       = zeros(Float64, NV2, NV2)         
const TAB_dS3_dy2       = zeros(Float64, NV2, NV2)  

##################################################
# Lyapunov NL part
##################################################

const TAB_DIV_v1_I1 = zeros(Float64, NV2)
const TAB_NL_S3 = zeros(Float64, NV2, NV2, NV2)

##################################################
# Lyapunov Source part
##################################################
const TAB_SOURCE_S3 = zeros(Float64, NV2, NV2, NV2)

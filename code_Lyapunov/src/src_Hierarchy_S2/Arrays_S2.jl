##################################################
# PREALLOCATED ARRAYS
##################################################

##################################################
# Lyapunov part
##################################################

# ATTENTION, here TAB_S2 correspond to S_2^(1) and not 1/N S_2^(1)
const TAB_S2            = zeros(Float64, NV2, NV2)    # 2 point correlation fonction 
const TAB_S2_SLICE      = zeros(Float64, NV2)         # S(v, v' fixed) 
const TAB_dS2_SLICE_dx  = zeros(Float64, NV2)         #      
const TAB_dS2_SLICE_dy  = zeros(Float64, NV2)         #
const TAB_dS2_dt        = zeros(Float64, NV2, NV2)   
#
const TAB_dS2_dx1        = zeros(Float64, NV2, NV2)              
const TAB_dS2_dy1        = zeros(Float64, NV2, NV2)   
const TAB_dS2_dx2       = zeros(Float64, NV2, NV2)         
const TAB_dS2_dy2       = zeros(Float64, NV2, NV2)  
#
const TAB_S2_SLICE_PARALLEL = [zeros(Float64, NV2) for tid = 1:Threads.nthreads()]
const TAB_DS2_DX_PARALLEL   = [zeros(Float64, NV2) for tid = 1:Threads.nthreads()]
const TAB_DS2_DY_PARALLEL   = [zeros(Float64, NV2) for tid = 1:Threads.nthreads()]
##################################################
# Lyapunov Source term 
##################################################

const TAB_Q = zeros(Float64, NV2, NV2) # Change the name 

##################################################
# Observables and Invariants
##################################################
const TAB_VAR_M   = zeros(Float64, NSTEPS+1) # Variance of total system mass

const TAB_VAR_Tx  = zeros(Float64, NSTEPS+1) # Variance of T_x(t)
const TAB_VAR_Ty  = zeros(Float64, NSTEPS+1) # Variance of T_y(t)
const TAB_VAR_E   = zeros(Float64, NSTEPS+1) # Variance of total system kinetic energy

const TAB_VAR_Px  = zeros(Float64, NSTEPS+1)   # Variance of total momentum x-direction
const TAB_VAR_Py  = zeros(Float64, NSTEPS+1)   # Variance of total momentum y-direction
const TAB_VAR_Ptot  = zeros(Float64, NSTEPS+1) # Variance of total momentum

const TAB_M4_Tx   = zeros(Float64, NSTEPS+1) 
const TAB_M4_Ty   = zeros(Float64, NSTEPS+1)  

# Cubic Observable VAR arrays
const TAB_VAR_SXXX = zeros(Float64,NSTEPS+1)
const TAB_VAR_SXXY = zeros(Float64,NSTEPS+1)
const TAB_VAR_SXYY = zeros(Float64,NSTEPS+1)
const TAB_VAR_SYYY = zeros(Float64,NSTEPS+1)
#
const TAB_VAR_NON_POLY_OBS = zeros(Float64,NSTEPS+1)

# Quartic Observable VAR arrays
const TAB_VAR_QXXXX = zeros(Float64,NSTEPS+1)
const TAB_VAR_QXXXY = zeros(Float64,NSTEPS+1)
const TAB_VAR_QXXYY = zeros(Float64,NSTEPS+1)
const TAB_VAR_QXYYY = zeros(Float64,NSTEPS+1)
const TAB_VAR_QYYYY = zeros(Float64,NSTEPS+1)
# Poly-sin-cos and Gaussian observables VAR arrays
const TAB_VAR_POLYSINCOS = zeros(Float64,NSTEPS+1)
const TAB_VAR_GAUSSIAN   = zeros(Float64,NSTEPS+1)
# Anisotropic cos channel pair, vsin_vsin, tanh/sigmoid observables VAR arrays
const TAB_VAR_COS_PLUS    = zeros(Float64,NSTEPS+1)
const TAB_VAR_COS_MINUS   = zeros(Float64,NSTEPS+1)
const TAB_VAR_VSIN_VSIN   = zeros(Float64,NSTEPS+1)
const TAB_VAR_TANH_VX     = zeros(Float64,NSTEPS+1)
const TAB_VAR_SIGMOID_VX  = zeros(Float64,NSTEPS+1)
const TAB_VAR_TANH_ANISO  = zeros(Float64,NSTEPS+1)

####################################################################################################
# Precompute integrals with F and wto S_2
# These arrays hold int dv3 B(v1, v3) F(v3) and int dv3 B(v1, v3) p_v3F(v3)
####################################################################################################
const TAB_INT_Bxx_F   = zeros(Float64, NV2)
const TAB_INT_Bxy_F   = zeros(Float64, NV2)
const TAB_INT_Byy_F   = zeros(Float64, NV2)
#
const TAB_INT_Bxx_dFx = zeros(Float64, NV2)
const TAB_INT_Bxy_dFx = zeros(Float64, NV2)
const TAB_INT_Bxy_dFy = zeros(Float64, NV2)
const TAB_INT_Byy_dFy = zeros(Float64, NV2)
#
const TAB_Ix_PARALLEL       = [zeros(Float64, NV2) for tid = 1:Threads.nthreads()]
const TAB_Iy_PARALLEL       = [zeros(Float64, NV2) for tid = 1:Threads.nthreads()]
const TAB_DIV_I_PARALLEL    = [zeros(Float64, NV2) for tid = 1:Threads.nthreads()]
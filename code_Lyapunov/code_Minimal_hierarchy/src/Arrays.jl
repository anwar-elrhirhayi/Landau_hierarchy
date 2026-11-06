##################################################
# PREALLOCATED ARRAYS
##################################################

##################################################
# VELOCITY GRID PARAMETERS
##################################################

# Volume element = DV^2
const DV = 2 * VMAX / NV  # DV_x = DV_y = DV
const NV2 = NV * NV
const DV2 = DV * DV
# Velocity grid
const TAB_Vx = TAB_Vx = [ -VMAX + (i - 0.5) * DV for i in 1:NV] # vx grid
const TAB_Vy = copy(TAB_Vx)                                     # vy grid

##################################################
# CORRESPONDENCE (i, j) --> l
##################################################

const TAB_i = [div(l - 1, NV) + 1 for l in 1:NV2]
const TAB_j = [(l - 1) % NV + 1 for l in 1:NV2]

##################################################
# Collision kernel
##################################################

const TAB_Bxx = zeros(Float64, NV2, NV2) 
const TAB_Bxy = zeros(Float64, NV2, NV2) 
const TAB_Byy = zeros(Float64, NV2, NV2) 

##################################################
# Observables and Invariants
##################################################

const TAB_TIME     = zeros(Float64,    NSTEPS+1) # Time array
const TAB_Tx       = zeros(Float64,    NSTEPS+1) # Temperature T_x(t)
const TAB_Ty       = zeros(Float64,    NSTEPS+1) # Temperature T_y(t)
const TAB_E        = zeros(Float64,    NSTEPS+1) # Total system kinetic energy
const TAB_Px       = zeros(Float64,    NSTEPS+1) # Total momentum x-direction
const TAB_Py       = zeros(Float64,    NSTEPS+1) # Total momentum y-direction
const TAB_M        = zeros(Float64,    NSTEPS+1) # Total system mass
const TAB_OBSERV   = zeros(Float64,    NSTEPS+1) # Observable used for the NV-convergence tests
# Cubic Observable AVG arrays
const TAB_SXXX = zeros(Float64,NSTEPS+1)
const TAB_SXXY = zeros(Float64,NSTEPS+1)
const TAB_SXYY = zeros(Float64,NSTEPS+1)
const TAB_SYYY = zeros(Float64,NSTEPS+1)
# Additionnal non polynomial observables
const TAB_NON_POLY_OBS = zeros(Float64,NSTEPS+1)
# Quartic Observable AVG arrays
const TAB_QXXXX = zeros(Float64,NSTEPS+1)
const TAB_QXXXY = zeros(Float64,NSTEPS+1)
const TAB_QXXYY = zeros(Float64,NSTEPS+1)
const TAB_QXYYY = zeros(Float64,NSTEPS+1)
const TAB_QYYYY = zeros(Float64,NSTEPS+1)
# Poly-sin-cos and Gaussian observables
const TAB_POLYSINCOS = zeros(Float64,NSTEPS+1)
const TAB_GAUSSIAN   = zeros(Float64,NSTEPS+1)
# Anisotropic cos channel pair, vsin_vsin, tanh/sigmoid observables
const TAB_COS_PLUS    = zeros(Float64,NSTEPS+1)
const TAB_COS_MINUS   = zeros(Float64,NSTEPS+1)
const TAB_VSIN_VSIN   = zeros(Float64,NSTEPS+1)
const TAB_TANH_VX     = zeros(Float64,NSTEPS+1)
const TAB_SIGMOID_VX  = zeros(Float64,NSTEPS+1)
const TAB_TANH_ANISO  = zeros(Float64,NSTEPS+1)

##################################################
# Landau eq part
##################################################

const TAB_F         = zeros(Float64, NV2)    # Distribution function F(v_x,v_y)
const TAB_dF_dx     = zeros(Float64, NV2)    # dF/dv_x
const TAB_dF_dy     = zeros(Float64, NV2)    # dF/dv_y
const TAB_dF_dt     = zeros(Float64, NV2)    # RHS of p_t F = div Flux
const TAB_Fluxx     = zeros(Float64, NV2)    # x-direction Flux 
const TAB_Fluxy     = zeros(Float64, NV2)    # y-direction Flux
const TAB_dFluxx_dx = zeros(Float64, NV2)    # dFluxx/dvx
const TAB_dFluxy_dy = zeros(Float64, NV2)    # dFluxy/dvy


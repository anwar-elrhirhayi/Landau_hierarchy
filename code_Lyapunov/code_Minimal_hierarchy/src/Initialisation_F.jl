##################################################
# Initial conditions for F
# Anisotropic Maxwellian distribution normalized to MASS_INIT
##################################################
function init_F_MAXW!(Tx, Ty)
    #### Fill TAB_F with anisotropic Maxwellian
    C_x = sqrt(1 / (2 * pi * Tx))
    C_y = sqrt(1 / (2 * pi * Ty))
    ####
    for l = 1:NV2
        i, j = TAB_i[l], TAB_j[l]
        vy = TAB_Vy[j]
        vx = TAB_Vx[i]
        TAB_F[l] = C_x * C_y * exp(-0.5 * (vx^2 / Tx + vy^2 / Ty))
    end
end

##################################################
# Initial conditions for F
# DOOR distribution normalized to MASS_INIT = 1
##################################################

# # 1D DF
# function gaussian_door(vx, a, n)    
#     norm = n / (a * gamma(1 / (2 * n)))
#     return norm * exp(-(vx/a)^(2 * n))
# end

# # F(t=0) smooth DOOR
# function init_F_DOOR_old!(Tx, Ty)
#     n = 8 # ATTENTION, for n > 9 the DOOR becomes very sharp. Use n = 8
#     alpha = 1/(2 * n)
#     a_x = sqrt(Tx * gamma(alpha) / gamma(3 * alpha))
#     a_y = sqrt(Ty * gamma(alpha) / gamma(3 * alpha))
        
#     for l = 1:NV2
#         #####
#         i, j = TAB_i[l], TAB_j[l]
#         vx = TAB_Vx[i]
#         vy = TAB_Vy[j]
#         #####
#         f_x = gaussian_door(vx, a_x, n)
#         f_y = gaussian_door(vy, a_y, n)
#         #####
#         TAB_F[l] = f_x * f_y
#         # ATTENTION -- F here is already normalized
#     end
# end

##################################################
# 1D smooth door obtained by convolving a uniform distribution on [-b, b]
# with a centred Gaussian
##################################################
function convolved_gaussian_door(v, b, sigma)
    # Sigma = Sharpness
    # b = width
    sqrt2_sigma = sqrt(2.0) * sigma
    return (erf((v + b) / sqrt2_sigma) - erf((v - b) / sqrt2_sigma)) / (4.0 * b)
    # ATTENTION -- thid distribution is already normalized
end

##################################################
# F(t=0): product of two Gaussian-convolved uniform distributions.
# For each direction, Var(U + Z) = b^2 / 3 + sigma^2 = T.
##################################################
function init_F_DOOR!(Tx, Ty)
    sigma = 0.15 # Sharpness parameter
    # ATTENTION -- for sigma < 0.14, the door becomes very sharp and
    # the initial conditions fluctuates a lot when NV changes
    # For sigma = 0.15, stabilisation at NV = 40
    if sigma^2 >= min(Tx, Ty)
        error("init_F_DOOR! -- sigma^2 must be smaller than Tx and Ty")
    end
    # Width b_x st Var v_x^2 = Tx at t=0
    # Idem for b_y
    b_x = sqrt(3.0 * (Tx - sigma^2))
    b_y = sqrt(3.0 * (Ty - sigma^2))

    @inbounds for l in 1:NV2
        #####
        i, j = TAB_i[l], TAB_j[l]
        vx = TAB_Vx[i]
        vy = TAB_Vy[j]

        f_x = convolved_gaussian_door(vx, b_x, sigma)
        f_y = convolved_gaussian_door(vy, b_y, sigma)
        TAB_F[l] = f_x * f_y
        # ATTENTION -- no normalisation is needed here
        #####
    end
end

##################################################
# Initial conditions for F
##################################################

if IC == "MAXW"
    const init_F! = init_F_MAXW!
elseif IC == "DOOR"
    const init_F! = init_F_DOOR!
end

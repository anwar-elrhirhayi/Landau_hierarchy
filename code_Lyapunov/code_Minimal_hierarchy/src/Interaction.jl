##################################################
# Definition of ln_lambda
# ATTENTION, this is for dimension 2
##################################################
function ln_lambda()
    return 1/(KMIN) - 1/(KMAX)
end
#####
const LN_LAMBDA = ln_lambda() # Computing once lnLambda
##################################################
# Function f in the definition of sigma
# It encodes the interaction potential
##################################################
function g_sigma(u)
    return 4 * pi * G^2 * LN_LAMBDA * MASS_PART / u^(3)
end
##################################################
# Kernel B(u) with u = v - v'
# B_{ij} = C * (|u|^2 Delta_ij - u_i u_j) / (|u|^2 + eps^2)^(3/2)
##################################################
function get_kernel_B(ux, uy)
    #####
    ux_SQ = ux * ux
    uy_SQ = uy * uy
    #####
    r2 = ux_SQ + uy_SQ # r^2 = ||v - v'||^2
    #####
    # We deal specifically with the exact resonance,
    # where the two velocities match
    if r2 == 0.0
        #####
        return 0.0, 0.0, 0.0, 0.0
        #####
    elseif r2 != 0
        #####
        r   = sqrt(r2)
        C   =   g_sigma(r)
        Bxx =   C * uy_SQ
        Byy =   C * ux_SQ
        Bxy = - C * ux * uy
        Byx = Bxy
        #####
        return Bxx, Bxy, Byx, Byy
    end
end
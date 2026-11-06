##################################################
# Covariance kernel C(v1, v2 )
# C(v1, v2) = p_v1 p_v2 . [ -F(v)F(v2 ) B(v,v2)
#            + delta(v1-v2) F(v1) int dv3 F(v3) B(v1,v3)]
##################################################


##################################################
# Compute the flux term in the covariance kernel Q:
# FLUX_Q(v1, v2) = - F(v1) F(v2) B(v1, v2) 
#                + delta(v1-v2) F(v1) INT_B_F(v1)
##################################################
const TAB_FLUX_Qxx = zeros(Float64, NV2, NV2)
const TAB_FLUX_Qxy = zeros(Float64, NV2, NV2)
const TAB_FLUX_Qyy = zeros(Float64, NV2, NV2)

function TAB_FLUX_Q!()
    #####
    fill!(TAB_FLUX_Qxx, 0.0)
    fill!(TAB_FLUX_Qxy, 0.0)
    fill!(TAB_FLUX_Qyy, 0.0)

    Threads.@threads for l1 in 1:NV2
        #####
        F_l1  =   TAB_F[l1]
        int_F_Bxx = TAB_INT_Bxx_F[l1]
        int_F_Bxy = TAB_INT_Bxy_F[l1]
        int_F_Byy = TAB_INT_Byy_F[l1]
        #####
        for l2 in 1:NV2
            #####
            F_l2 = TAB_F[l2]
            #####
            # Read B(v1, v2)
            Bxx = TAB_Bxx[l1, l2]
            Bxy = TAB_Bxy[l1, l2]
            Byy = TAB_Byy[l1, l2]
            ####
            delta = delta_cell(l1, l2)
            TAB_FLUX_Qxx[l1, l2] = -F_l1 * F_l2 * Bxx + delta * F_l1 * int_F_Bxx
            TAB_FLUX_Qxy[l1, l2] = -F_l1 * F_l2 * Bxy + delta * F_l1 * int_F_Bxy
            TAB_FLUX_Qyy[l1, l2] = -F_l1 * F_l2 * Byy + delta * F_l1 * int_F_Byy
            ####
        end
    end
end

##################################################
# Compute p_v2 . FLUX_Q(v1,v2)
# by finite differences on indices l2
##################################################
const TAB_DIV_FLUX_Qx = zeros(Float64, NV2, NV2)
const TAB_DIV_FLUX_Qy = zeros(Float64, NV2, NV2)

function div_prime_FLUX_Q!()
    fill!(TAB_DIV_FLUX_Qx, 0.0)
    fill!(TAB_DIV_FLUX_Qy, 0.0)
    ####
    # Div_FLUX_Q = (p_v2x, p_v2y)^T * FLUX_Q 
    ####
    Threads.@threads for l1 in 1:NV2
        for l2 in 1:NV2
            ####
            # DIV_FLUX_Qx = p_v2x FLUX_Qxx + p_v2y FLUX_Qyx
            ####
            dFLUX_Qxx = value_div_2D_field_v2(TAB_FLUX_Qxx, 1, l1, l2)
            dFLUX_Qxy = value_div_2D_field_v2(TAB_FLUX_Qxy, 2, l1, l2)

            TAB_DIV_FLUX_Qx[l1, l2] = dFLUX_Qxx + dFLUX_Qxy

            ####
            # DIV_FLUX_Qy = p_v2x FLUX_Qxy + p_v2y FLUX_Qyy
            ####
            dFLUX_Qyx = value_div_2D_field_v2(TAB_FLUX_Qxy, 1, l1, l2)
            dFLUX_Qyy = value_div_2D_field_v2(TAB_FLUX_Qyy, 2, l1, l2)
            # By symmetry of B, dFLUX_Qxy = dFLUX_Qyx
            TAB_DIV_FLUX_Qy[l1, l2] = dFLUX_Qyx + dFLUX_Qyy
            ####
        end
    end
end


##################################################
# Now apply p_v1 to the result to get C
##################################################
function compute_Source_S2!()
    fill!(TAB_SOURCE_S2, 0.0)
    #####
    TAB_FLUX_Q!()       # Compute the flux (term inside \p_v1 \p_v2 )
    div_prime_FLUX_Q!() # Compute \p_v2  flux_Q
    #####
    for l1 in 1:NV2
        for l2 in 1:l1
            ####
            # p_v1x DIV_FLUX_Qx
            ####
            ddFLUX_Qxx2 = value_div_2D_field_v1(TAB_DIV_FLUX_Qx, 1, l1, l2)
            ####
            # p_vy DIV_FLUX_Qy 
            ####
            ddFLUX_Qyy2 = value_div_2D_field_v1(TAB_DIV_FLUX_Qy, 2, l1, l2)
            ####
            # p_v1x [p_v2x FL_Qxx + p_v2y FL_Qyx] 
            #  + p_v1y [p_v2x FL_Qxy + p_v2y FL_Qyy] 
            TAB_SOURCE_S2[l1, l2] += (ddFLUX_Qxx2 + ddFLUX_Qyy2)
            TAB_SOURCE_S2[l2, l1] = TAB_SOURCE_S2[l1, l2]
        end
    end
end

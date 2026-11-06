
##################################################
# Compute (v1, v2) --> int dv S2(v1, v) B(v, v2)
##################################################
const TAB_INT_S2_Bxx = zeros(Float64, NV2, NV2)
const TAB_INT_S2_Bxy = zeros(Float64, NV2, NV2)
const TAB_INT_S2_Byy = zeros(Float64, NV2, NV2)

function precompute_S2_B_integral!()
    fill!(TAB_INT_S2_Bxx, 0.0)
    fill!(TAB_INT_S2_Bxy, 0.0)
    fill!(TAB_INT_S2_Byy, 0.0)
    #####
    Threads.@threads for l1 in 1:NV2
        for l2 in 1:NV2
            #####
            acc11 = 0.0
            acc12 = 0.0
            acc22 = 0.0
            #####
            for l in 1:NV2
                #####
                S2_l1_l = TAB_S2[l1, l]
                # Read the component in TAB_B
                B11 = TAB_Bxx[l, l2]
                B12 = TAB_Bxy[l, l2]
                B22 = TAB_Byy[l, l2]
                #####
                acc11 += S2_l1_l * B11
                acc12 += S2_l1_l * B12
                acc22 += S2_l1_l * B22
                #####
            end
            #####
            TAB_INT_S2_Bxx[l1, l2] = acc11 * DV2
            TAB_INT_S2_Bxy[l1, l2] = acc12 * DV2
            TAB_INT_S2_Byy[l1, l2] = acc22 * DV2
        end
    end
end

####################################################################################################
# Precompute integrals with F in the Linear part of Lyapunov for S2 and S3
# These arrays hold int dv B(v1, v) F(v) and int dv B(v1, v) p_v F(v)
####################################################################################################


##################################################
# Precompute the integrals for Linear part of Lyapunov
# v1 --> TAB_INT_B_F(v1)  = int dv B(v1, v) F(v)
# v1 --> TAB_INT_B_dF(v1) = int dv B(v1, v) p_v F(v)
##################################################
function precompute_F_B_integrals!()
    # Loop over v1 = l1
    Threads.@threads :static for l1 in 1:NV2
        #####
        # Read the collision kernel slices at fixed v1 ie B(v1, .)
        tab_Bxx = TAB_TAB_Bxx[l1]
        tab_Bxy = TAB_TAB_Bxy[l1]
        tab_Byy = TAB_TAB_Byy[l1]
        #
        acc_Bxx_F   = 0.0
        acc_Bxy_F   = 0.0
        acc_Byy_F   = 0.0
        #
        acc_Bxx_dFx = 0.0
        acc_Bxy_dFx = 0.0
        acc_Bxy_dFy = 0.0
        acc_Byy_dFy = 0.0

        #####
        # Loop over v = l
        #####
        @inbounds @simd for l in 1:NV2
            # Values of F and its gradients at v
            F     = TAB_F[l]
            dF_dx = TAB_dF_dx[l]
            dF_dy = TAB_dF_dy[l]
            # Values of B(v1, v)
            Bxx = tab_Bxx[l]
            Bxy = tab_Bxy[l]
            Byy = tab_Byy[l]

            # Compute the integral int dv B(v1, v) F(v)
            acc_Bxx_F   += @muladd(Bxx * F)
            acc_Bxy_F   += @muladd(Bxy * F)
            # Here we use the symmetry acc_B_yx_F = acc_Bxy_F
            acc_Byy_F   += @muladd(Byy * F)

            # Compute the integral int dv3 B(v1, v) p_v F(v)
            acc_Bxx_dFx += @muladd(Bxx * dF_dx)
            acc_Bxy_dFx += @muladd(Bxy * dF_dx)
            acc_Bxy_dFy += @muladd(Bxy * dF_dy)
            acc_Byy_dFy += @muladd(Byy * dF_dy)
        end
        
        #####
        # Filling in the global arrays
        TAB_INT_Bxx_F[l1]   = acc_Bxx_F * DV2
        TAB_INT_Bxy_F[l1]   = acc_Bxy_F * DV2
        TAB_INT_Byy_F[l1]   = acc_Byy_F * DV2
        #
        TAB_INT_Bxx_dFx[l1] = acc_Bxx_dFx * DV2
        TAB_INT_Bxy_dFx[l1] = acc_Bxy_dFx * DV2
        TAB_INT_Bxy_dFy[l1] = acc_Bxy_dFy * DV2
        TAB_INT_Byy_dFy[l1] = acc_Byy_dFy * DV2
    end
end

##################################################
##################################################
# Integrals arrays
##################################################

const TAB_INT_I2_xx = zeros(Float64, NV2, NV2)
const TAB_INT_I2_xy = zeros(Float64, NV2, NV2)
const TAB_INT_I2_yy = zeros(Float64, NV2, NV2)

const TAB_INT_I3_x = zeros(Float64, NV2, NV2)
const TAB_INT_I3_y = zeros(Float64, NV2, NV2)

##################################################
# Compute the 3 integral blocks for the NL term
##################################################
function TAB_INTEGRALS_NL!()
    
    fill!(TAB_INT_I2_xx, 0.0)
    fill!(TAB_INT_I2_xy, 0.0)
    fill!(TAB_INT_I2_yy, 0.0)
    
    fill!(TAB_INT_I3_x, 0.0)
    fill!(TAB_INT_I3_y, 0.0)
    
    # Computation for I2 (Tensor) and I3 (Vector)
    Threads.@threads for l1 in 1:NV2
        for l2 in 1:NV2
            ####
            acc_I2_xx = 0.0
            acc_I2_xy = 0.0
            acc_I2_yy = 0.0
            #
            acc_I3_x  = 0.0
            acc_I3_y  = 0.0
            ####
            for l in 1:NV2
                Bxx = TAB_Bxx[l1, l]
                Bxy = TAB_Bxy[l1, l]
                Byy = TAB_Byy[l1, l]
                
                # I2 term: int dv B(v1, v) S2(v, v2)
                S2_v_v2 = TAB_S2[l, l2]
                acc_I2_xx += Bxx * S2_v_v2
                acc_I2_xy += Bxy * S2_v_v2
                acc_I2_yy += Byy * S2_v_v2
                
                # I3 term: int dv B(v1, v) \p_v S2(v, v2)
                gradS2_x = TAB_dS2_dx[l, l2]
                gradS2_y = TAB_dS2_dy[l, l2]
                acc_I3_x += (Bxx * gradS2_x + Bxy * gradS2_y)
                acc_I3_y += (Bxy * gradS2_x + Byy * gradS2_y)
            end
            
            TAB_INT_I2_xx[l1, l2] = acc_I2_xx * DV2
            TAB_INT_I2_xy[l1, l2] = acc_I2_xy * DV2
            TAB_INT_I2_yy[l1, l2] = acc_I2_yy * DV2
            
            TAB_INT_I3_x[l1, l2] = acc_I3_x * DV2
            TAB_INT_I3_y[l1, l2] = acc_I3_y * DV2
        end
    end
end
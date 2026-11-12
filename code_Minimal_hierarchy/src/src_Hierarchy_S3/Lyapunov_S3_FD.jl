####################################################################################################
# PDE for S(v1, v2,t) for Landau - finite differences
####################################################################################################

####################################################################################################
# For a given TAB_S2, compute contribution L_v1 S(v1, v2, v3) if indx == 0
#                                          L_v2 S(v1, v2, v3) if indx == 1
#                                          L_v3 S(v1, v2, v3) if indx == 2
####################################################################################################


##################################################
# Compute Lv
##################################################

function compute_Lv_contribution_S3!(TAB_F, TAB_dF_dx, TAB_dF_dy, indx)
    # Loop over v2 and v3 --> index l2, l3
    Threads.@threads :static  for l3 = 1:NV2
        for l2 in 1:NV2 
            #####
            tid = Threads.threadid()
            #####
            # Allocate thread-local temporary arrays to avoid data races
            # (each thread has its own copies of these buffers)
            # without new allocations
            tab_S3_slice     = TAB_S3_SLICE_PARALLEL[tid]
            tab_dS3_dx_slice = TAB_DS3_DX_PARALLEL[tid]
            tab_dS3_dy_slice = TAB_DS3_DY_PARALLEL[tid]
            Integrand_x      = TAB_Ix_PARALLEL[tid]
            Integrand_y      = TAB_Iy_PARALLEL[tid]
            div_Int          = TAB_DIV_I_PARALLEL[tid]

            #####
            # From here, v2 = l2, v3 = l3 is fixed.
            # Build S3_slice(v1) = S3(v1, v2, v3) if indx == 0 or
            #                    = S3(v2, v1, v3) if indx == 1 or
            #                    = S3(v2, v3, v1) if indx == 2 or

            # Each thread calls this with its own buffers.
            slice_3D_field(indx, l2, l3, tab_S3_slice)

            #####
            # Compute gradients with respect to v1 for S2_slice
            # gradient = p_v1 S3(v1, v2, v3) if indx == 0 or
            #          = p_v1 S3(v2, v1, v3) if indx == 1 or
            #          = p_v1 S3(v2, v3, v1) if indx == 2

            # dS2_dx_slice = p_v1x S3_slice(v1)
            # dS2_dy_slice = p_v1y S3_slice(v1)
            compute_1D_field_gradients!(tab_S3_slice, tab_dS3_dx_slice, tab_dS3_dy_slice)
            
            #####
            # For each v1 = l1, accumulate the integral over v = l
            # All operations here are thread-local, so no synchronization needed
            @inbounds for l1 in 1:NV2
                # Values of F and its gradients at v1
                Fv    = TAB_F[l1]
                dF_dx = TAB_dF_dx[l1]
                dF_dy = TAB_dF_dy[l1]
                
                # Local S3(v1, v2, v3)
                S3_v  = tab_S3_slice[l1] 
                # p_v1 S3(v1, v2, v3) for v2, v3 fixed
                dS3_x = tab_dS3_dx_slice[l1]
                dS3_y = tab_dS3_dy_slice[l1]
                
                # Collision kernel B(v, v1)
                # Tab_Bab = B_ab(v1, .) at fixed v1
                tab_Bxx = TAB_TAB_Bxx[l1]
                tab_Bxy = TAB_TAB_Bxy[l1] 
                tab_Byy = TAB_TAB_Byy[l1]
                
                # Accumulators for the integrals that depends on S_3
                int_Bxx_S3   = 0.0
                int_Bxy_S3   = 0.0
                int_Byy_S3   = 0.0
                #
                int_Bxx_dS3x = 0.0
                int_Bxy_dS3x = 0.0
                int_Bxy_dS3y = 0.0
                int_Byy_dS3y = 0.0

                #####
                # For each v1, v2, v3, we accumulate the integral over v = l
                #####
                @inbounds @simd for l in 1:NV2
                    # S_3 at v (for fixed v2, v3) = S3(v, v2, v3) 
                    S3_l   = tab_S3_slice[l]
                    # Gradient of S3 at v ie p_v S(v, v2, v3)
                    dS3x_l = tab_dS3_dx_slice[l]
                    dS3y_l = tab_dS3_dy_slice[l]

                    # Collision kernel B(v, v1)
                    Bxx = tab_Bxx[l]
                    Bxy = tab_Bxy[l] 
                    Byy = tab_Byy[l]

                    # Compute the integral int dv B(v, v1) S3(v, v2, v3) 
                    int_Bxx_S3   = @muladd(int_Bxx_S3 + Bxx * S3_l)
                    int_Bxy_S3   = @muladd(int_Bxy_S3 + Bxy * S3_l)
                    int_Byy_S3   = @muladd(int_Byy_S3 + Byy * S3_l)
                    # Compute the integral int dv B(v1, v) p_v S2(v, v2, v3) 
                    int_Bxx_dS3x = @muladd(int_Bxx_dS3x + Bxx * dS3x_l)
                    int_Bxy_dS3x = @muladd(int_Bxy_dS3x + Bxy * dS3x_l)
                    int_Bxy_dS3y = @muladd(int_Bxy_dS3y + Bxy * dS3y_l)
                    int_Byy_dS3y = @muladd(int_Byy_dS3y + Byy * dS3y_l)
                end

                #####
                # Groups all the terms
                # term = p_v1 F(v1) int_B_S3 - S_3(v1, v2, v3) int_B_p_v1F + p_v S3(v1,v2,v3) int_B_F - F(v) int_B_p_v1 S3
                #####
                term_xx = @muladd((dF_dx * int_Bxx_S3 * DV2) - S3_v * TAB_INT_Bxx_dFx[l1] + dS3_x * TAB_INT_Bxx_F[l1] - (Fv * int_Bxx_dS3x * DV2))
                term_xy = @muladd((dF_dy * int_Bxy_S3 * DV2) - S3_v * TAB_INT_Bxy_dFy[l1] + dS3_y * TAB_INT_Bxy_F[l1] - (Fv * int_Bxy_dS3y * DV2))
                acc_x = term_xx + term_xy
                #
                term_yx = @muladd((dF_dx * int_Bxy_S3 * DV2) - S3_v * TAB_INT_Bxy_dFx[l1] + dS3_x * TAB_INT_Bxy_F[l1] - (Fv * int_Bxy_dS3x * DV2))
                term_yy = @muladd((dF_dy * int_Byy_S3 * DV2) - S3_v * TAB_INT_Byy_dFy[l1] + dS3_y * TAB_INT_Byy_F[l1] - (Fv * int_Byy_dS3y * DV2))
                acc_y = term_yx + term_yy

                # Accumulators for the integral result
                # ATTENTION to the factor 1/2
                Integrand_x[l1] = 0.5 * acc_x 
                Integrand_y[l1] = 0.5 * acc_y
            end

            #####
            # Compute the divergence of TAB_Integrand (p_v1x Integrand_x + p_v1y Integrand_y)
            # by finite difference.
            # Gives p_tS contribution for this v2 = l2, v3 = l3
            compute_1D_field_div!(Integrand_x, Integrand_y, div_Int)
            #####
            # Accumulate this (v2, v3) contribution into the global array.
            # ATTENTION, the accumulation depends on the direction = v1, v2 or v3
            if indx == 0     # We are computing L_v1 S_3(v1, v2, v3)
                for l1 in 1:NV2
                    TAB_dS3_dt[l1, l2, l3] += div_Int[l1]
                end
            #
            elseif indx == 1 # We are computing L_v2 S_3(v3, v1, v2)
                for l1 in 1:NV2
                    TAB_dS3_dt[l3, l1, l2] += div_Int[l1]
                end
            #    
            elseif indx == 2 # We are computing L_v3 S_3(v2, v3, v1)
                for l1 in 1:NV2
                    TAB_dS3_dt[l2, l3, l1] += div_Int[l1]
                end
            end
        end
    end
end


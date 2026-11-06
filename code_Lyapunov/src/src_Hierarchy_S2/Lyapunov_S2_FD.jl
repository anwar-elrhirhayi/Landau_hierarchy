####################################################################################################
# PDE for S(v1, v2,t) for Landau - finite differences
####################################################################################################

####################################################################################################
# For a given TAB_S2, compute contribution L_v  S(v1, v2) if indx == 0
#                                          L_v2 S(v1, v2) if indx == 1
# In practice,
# If indx == 0, S2_SLICE(v) = S(v2, v1)
# If indx == 1, S2_SLICE(v) = S(v1, v2)
# and then we compute L_v1 S2_SLICE(v).
# More precisely, we compute L_v1 S2_SLICE = p_v1 I(v1).
#   I_x(l1) = int dv3 B(v1, v3) . [ ... ] (vx comp)
#   I_y(l1) = similar (y comp)
# then divergence of I:  L_v1 S2_SLICE = p_v1 I(v1).
# This function is then added into p_t S2.
####################################################################################################


####################################################################################################
# Lv - PARALLELISATION ON THE FIRST LOOP v2 = l2
####################################################################################################

##################################################
# Precompute the integrals
# TAB_INT_B_F = int dv3 B(v1, v3) F(v3)
# TAB_INT_B_dF = int dv3 B(v1, v3) p_v3 F(v3)
##################################################
function precompute_F_integrals!()
    # Loop over v1 = l1
    Threads.@threads :static for l1 in 1:NV2
        #####
        # Read the collision kernel slices at fixed v1 ie B(v1, .)
        tab_Bxx = TAB_TAB_Bxx[l1]
        tab_Bxy = TAB_TAB_Bxy[l1]
        tab_Byy = TAB_TAB_Byy[l1]

        acc_Bxx_F   = 0.0
        acc_Bxy_F   = 0.0
        acc_Byy_F   = 0.0
        #
        acc_Bxx_dFx = 0.0
        acc_Bxy_dFx = 0.0
        acc_Bxy_dFy = 0.0
        acc_Byy_dFy = 0.0

        #####
        # Loop over v3 = l3
        #####
        @inbounds @simd for l3 in 1:NV2
            # Values of F and its gradients at v3
            F3     = TAB_F[l3]
            dF_dx3 = TAB_dF_dx[l3]
            dF_dy3 = TAB_dF_dy[l3]
            # Values of B(v1, v3)
            Bxx = tab_Bxx[l3]
            Bxy = tab_Bxy[l3]
            Byy = tab_Byy[l3]

            # Compute the integral int dv3 B(v1, v3) F(v3)
            acc_Bxx_F   += @muladd(Bxx * F3)
            acc_Bxy_F   += @muladd(Bxy * F3)
            # Here we use the symmetry acc_B_yx_F = acc_Bxy_F
            acc_Byy_F   += @muladd(Byy * F3)

            # Compute the integral int dv3 B(v1, v3) p_v3 F(v3)
            acc_Bxx_dFx += @muladd(Bxx * dF_dx3)
            acc_Bxy_dFx += @muladd(Bxy * dF_dx3)
            acc_Bxy_dFy += @muladd(Bxy * dF_dy3)
            acc_Byy_dFy += @muladd(Byy * dF_dy3)
        end
        
        #####
        # Filling in the global arrays
        TAB_INT_Bxx_F[l1]   = acc_Bxx_F
        TAB_INT_Bxy_F[l1]   = acc_Bxy_F
        TAB_INT_Byy_F[l1]   = acc_Byy_F
        #
        TAB_INT_Bxx_dFx[l1] = acc_Bxx_dFx
        TAB_INT_Bxy_dFx[l1] = acc_Bxy_dFx
        TAB_INT_Bxy_dFy[l1] = acc_Bxy_dFy
        TAB_INT_Byy_dFy[l1] = acc_Byy_dFy
    end
end

####################################################################################################
# Compute L_v1 and L_v2
####################################################################################################

##################################################
# Compute_Lv using
# + Threads.@threads :static for l2
# + WITHOUT allocation
# + without the parallelisation on l1
##################################################
function compute_Lv_contribution!(TAB_F, TAB_dF_dx, TAB_dF_dy, indx)

    # Loop over v2 index pair (ip,jp)
    Threads.@threads :static for l2 in 1:NV2
        #####
        tid = Threads.threadid()
        #####
        # Allocate thread-local temporary arrays to avoid data races
        # (each thread has its own copies of these buffers)
        # without new allocations
        tab_S2_slice     = TAB_S2_SLICE_PARALLEL[tid]
        tab_dS2_dx_slice = TAB_DS2_DX_PARALLEL[tid]
        tab_dS2_dy_slice = TAB_DS2_DY_PARALLEL[tid]
        Integrand_x      = TAB_Ix_PARALLEL[tid]
        Integrand_y      = TAB_Iy_PARALLEL[tid]
        div_Int          = TAB_DIV_I_PARALLEL[tid]

        #####
        # From here, v2 = (i2, j2) is fixed.
        # Build S2_slice(v1) = S2(v1, v2) if indx == 0 or
        #                    = S2(v2, v1) if indx == 1.
        # Each thread calls this with its own buffers.
        slice_2D_field(indx, l2, tab_S2_slice)

        #####
        # Compute gradients with respect to v1 for S2_slice
        # gradient = p_v S2(v1, v2) if indx == 0 or
        #          = p_v S2(v2, v1) if indx == 1
        # dS2_dx_slice = p_v1x S2_slice(v1)
        # dS2_dy_slice = p_v1y S2_slice(v1)
        compute_1D_field_gradients!(tab_S2_slice, tab_dS2_dx_slice, tab_dS2_dy_slice)
        
        #####
        # For each v1 = l1, accumulate the integral over v3 = (ipp, jpp)
        # All operations here are thread-local1, so no synchronization needed
        @inbounds for l1 in 1:NV2
            # Values of F and its gradients at v
            Fv    = TAB_F[l1]
            dF_dx = TAB_dF_dx[l1]
            dF_dy = TAB_dF_dy[l1]
            
            # Local S2(v1, v2)
            S2_v  = tab_S2_slice[l1] 
            # p_v S2(v1, v2) for v2 fixed
            dS2_x = tab_dS2_dx_slice[l1]
            dS2_y = tab_dS2_dy_slice[l1]
            
            # Collision kernel B(v1, v3)
            # ATTENTION, B(v1, v3) and not B(v2, v3)
            # Tab_Bab = B_ab(v1, .) at fixed v1
            tab_Bxx = TAB_TAB_Bxx[l1]
            tab_Bxy = TAB_TAB_Bxy[l1] 
            tab_Byy = TAB_TAB_Byy[l1]
            
            # Accumulators for the integrals that depends on S_2
            int_Bxx_S2   = 0.0
            int_Bxy_S2   = 0.0
            int_Byy_S2   = 0.0
            #
            int_Bxx_dS2x = 0.0
            int_Bxy_dS2x = 0.0
            int_Bxy_dS2y = 0.0
            int_Byy_dS2y = 0.0

            #####
            # For each v1 = l1, accumulate the integral over v3 = l3
            #####
            @inbounds @simd for l3 in 1:NV2
                # S_2 at v3 (for fixed v2) = S2(v3, v2) 
                S2_l3   = tab_S2_slice[l3]
                # Gradient of S2 at v3 ie p_v3 S(v3, v2)
                dS2x_l3 = tab_dS2_dx_slice[l3]
                dS2y_l3 = tab_dS2_dy_slice[l3]

                # Collision kernel B(v1, v3)
                Bxx = tab_Bxx[l3]
                Bxy = tab_Bxy[l3] 
                Byy = tab_Byy[l3]

                # Compute the integral int dv3 B(v1, v3) S2(v3, v2) 
                int_Bxx_S2   = @muladd(int_Bxx_S2 + Bxx * S2_l3)
                int_Bxy_S2   = @muladd(int_Bxy_S2 + Bxy * S2_l3)
                int_Byy_S2   = @muladd(int_Byy_S2 + Byy * S2_l3)
                # Compute the integral int dv3 B(v1, v3) p_v3 S2(v3, v2) 
                int_Bxx_dS2x = @muladd(int_Bxx_dS2x + Bxx * dS2x_l3)
                int_Bxy_dS2x = @muladd(int_Bxy_dS2x + Bxy * dS2x_l3)
                int_Bxy_dS2y = @muladd(int_Bxy_dS2y + Bxy * dS2y_l3)
                int_Byy_dS2y = @muladd(int_Byy_dS2y + Byy * dS2y_l3)
            end

            #####
            # Groups all the terms
            # term = p_v1 F(v1) int_B_S2 - S_2(v1, v2) int_B_p_v1F + p_v S2(v1,v2) int_B_F - F(v) int_B_p_v1 S2
            #####
            term_xx = @muladd(dF_dx * int_Bxx_S2 - S2_v * TAB_INT_Bxx_dFx[l1] + dS2_x * TAB_INT_Bxx_F[l1] - Fv * int_Bxx_dS2x)
            term_xy = @muladd(dF_dy * int_Bxy_S2 - S2_v * TAB_INT_Bxy_dFy[l1] + dS2_y * TAB_INT_Bxy_F[l1] - Fv * int_Bxy_dS2y)
            acc_x = term_xx + term_xy
            #
            term_yx = @muladd(dF_dx * int_Bxy_S2 - S2_v * TAB_INT_Bxy_dFx[l1] + dS2_x * TAB_INT_Bxy_F[l1] - Fv * int_Bxy_dS2x)
            term_yy = @muladd(dF_dy * int_Byy_S2 - S2_v * TAB_INT_Byy_dFy[l1] + dS2_y * TAB_INT_Byy_F[l1] - Fv * int_Byy_dS2y)
            acc_y = term_yx + term_yy

            # Accumulators for the integral result
            Integrand_x[l1] = acc_x * DV2
            Integrand_y[l1] = acc_y * DV2
        end

        #####
        # Compute the divergence of TAB_Integrand (p_v1x Integrand_x + p_v1y Integrand_y)
        # by finite difference.
        # Gives p_tS contribution for this v2 = l2
        compute_1D_field_div!(Integrand_x, Integrand_y, div_Int)
        #####
        # Accumulate this v2 contribution into the global array.
        # ATTENTION, the accumulation depends on the direction = v1 or v2
        if indx == 0     # We are computing L_v S_2(v1, v2)
            for l1 in 1:NV2
                TAB_dS2_dt[l1, l2] += div_Int[l1]
            end
            #
        elseif indx == 1 # We are computing L_v2 S_2(v1, v2)
            for l1 in 1:NV2
                TAB_dS2_dt[l2, l1] += div_Int[l1]
            end
        end
    end
end


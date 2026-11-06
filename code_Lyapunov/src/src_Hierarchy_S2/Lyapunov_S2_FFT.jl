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
# Compute_Lv using
# + WITHOUT allocation
##################################################

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
# Containers for the parallelisation in L_v S2
##################################################
struct container_FFT_S2
    # Arrays in the R-space
    tab_S2_R::Array{Float64, 2}
    tab_dS2dx_R::Array{Float64, 2}
    tab_dS2dy_R::Array{Float64, 2}
    #
    tab_conv_Bxx_S2_R::Array{Float64, 2}
    tab_conv_Bxy_S2_R::Array{Float64, 2}
    tab_conv_Byy_S2_R::Array{Float64, 2}
    #
    tab_conv_Bxx_dS2dx_R::Array{Float64, 2}
    tab_conv_Bxy_dS2dx_R::Array{Float64, 2}
    tab_conv_Bxy_dS2dy_R::Array{Float64, 2}
    tab_conv_Byy_dS2dy_R::Array{Float64, 2}
    
    # Arrays in the K-space
    tab_S2_K::Array{Complex{Float64}, 2}
    tab_dS2dx_K::Array{Complex{Float64}, 2}
    tab_dS2dy_K::Array{Complex{Float64}, 2}
    #
    tab_conv_Bxx_S2_K::Array{Complex{Float64}, 2}
    tab_conv_Bxy_S2_K::Array{Complex{Float64}, 2}
    tab_conv_Byy_S2_K::Array{Complex{Float64}, 2}
    #
    tab_conv_Bxx_dS2dx_K::Array{Complex{Float64}, 2}
    tab_conv_Bxy_dS2dx_K::Array{Complex{Float64}, 2}
    tab_conv_Bxy_dS2dy_K::Array{Complex{Float64}, 2}
    tab_conv_Byy_dS2dy_K::Array{Complex{Float64}, 2}
end

const CONTAINER_FFT_S2 = [ container_FFT_S2(zeros(Float64,NV_DOUBLE,NV_DOUBLE),
                                            zeros(Float64,NV_DOUBLE,NV_DOUBLE),
                                            zeros(Float64,NV_DOUBLE,NV_DOUBLE),
                                            zeros(Float64,NV_DOUBLE,NV_DOUBLE),
                                            zeros(Float64,NV_DOUBLE,NV_DOUBLE),
                                            zeros(Float64,NV_DOUBLE,NV_DOUBLE),
                                            zeros(Float64,NV_DOUBLE,NV_DOUBLE),
                                            zeros(Float64,NV_DOUBLE,NV_DOUBLE),
                                            zeros(Float64,NV_DOUBLE,NV_DOUBLE),
                                            zeros(Float64,NV_DOUBLE,NV_DOUBLE),
                                            zeros(Complex{Float64},NV+1,NV_DOUBLE),
                                            zeros(Complex{Float64},NV+1,NV_DOUBLE),
                                            zeros(Complex{Float64},NV+1,NV_DOUBLE),
                                            zeros(Complex{Float64},NV+1,NV_DOUBLE),
                                            zeros(Complex{Float64},NV+1,NV_DOUBLE),
                                            zeros(Complex{Float64},NV+1,NV_DOUBLE),
                                            zeros(Complex{Float64},NV+1,NV_DOUBLE),
                                            zeros(Complex{Float64},NV+1,NV_DOUBLE),
                                            zeros(Complex{Float64},NV+1,NV_DOUBLE),
                                            zeros(Complex{Float64},NV+1,NV_DOUBLE)
                                            )
                                            for ithread=1:Threads.nthreads()]

##################################################
# Prepare the FFT of S2, p_v S2 
##################################################
function prepare_fields_S2!(tab_S2_slice, tab_dS2_dx_slice, tab_dS2_dy_slice, cnt_FFT_S2)
    #
    tab_S2_R     = cnt_FFT_S2.tab_S2_R
    tab_dS2dx_R  = cnt_FFT_S2.tab_dS2dx_R
    tab_dS2dy_R  = cnt_FFT_S2.tab_dS2dy_R
    #
    fill!(tab_S2_R   , 0.0)
    fill!(tab_dS2dx_R, 0.0)
    fill!(tab_dS2dy_R, 0.0)
    #
    @inbounds for l = 1:NV2
        #####
        i, j = TAB_i[l], TAB_j[l]
        #
        tab_S2_R[i,j]    = tab_S2_slice[l]
        tab_dS2dx_R[i,j] = tab_dS2_dx_slice[l]
        tab_dS2dy_R[i,j] = tab_dS2_dy_slice[l]
        #####
    end
end

##################################################
# Compute the FFT of S2 slices and its gradients
##################################################
function compute_FFT_S2!(cnt_FFT_S2)
    tab_S2_R        = cnt_FFT_S2.tab_S2_R
    tab_dS2dx_R     = cnt_FFT_S2.tab_dS2dx_R
    tab_dS2dy_R     = cnt_FFT_S2.tab_dS2dy_R
    #
    tab_S2_K        = cnt_FFT_S2.tab_S2_K
    tab_dS2dx_K     = cnt_FFT_S2.tab_dS2dx_K
    tab_dS2dy_K     = cnt_FFT_S2.tab_dS2dy_K

    # Compute FFT
    mul!(tab_S2_K,    PLAN_R_TO_K, tab_S2_R)
    mul!(tab_dS2dx_K, PLAN_R_TO_K, tab_dS2dx_R)
    mul!(tab_dS2dy_K, PLAN_R_TO_K, tab_dS2dy_R)
end

##################################################
# Compute the convolution in L_v S2
##################################################
function compute_convolution_S2!(cnt_FFT_S2)
    tab_S2_K        = cnt_FFT_S2.tab_S2_K
    tab_dS2dx_K     = cnt_FFT_S2.tab_dS2dx_K
    tab_dS2dy_K     = cnt_FFT_S2.tab_dS2dy_K
    #
    tab_conv_Bxx_S2_K = cnt_FFT_S2.tab_conv_Bxx_S2_K
    tab_conv_Bxy_S2_K = cnt_FFT_S2.tab_conv_Bxy_S2_K
    tab_conv_Byy_S2_K = cnt_FFT_S2.tab_conv_Byy_S2_K
    #
    tab_conv_Bxx_dS2dx_K = cnt_FFT_S2.tab_conv_Bxx_dS2dx_K
    tab_conv_Bxy_dS2dx_K = cnt_FFT_S2.tab_conv_Bxy_dS2dx_K
    tab_conv_Bxy_dS2dy_K = cnt_FFT_S2.tab_conv_Bxy_dS2dy_K
    tab_conv_Byy_dS2dy_K = cnt_FFT_S2.tab_conv_Byy_dS2dy_K

    # Compute the convolution
    @inbounds for j = 1:NV_DOUBLE
        @simd for i = 1:NV+1
            ####
            tab_conv_Bxx_S2_K[i,j] = TAB_Bxx_K[i,j] * tab_S2_K[i,j]
            tab_conv_Bxy_S2_K[i,j] = TAB_Bxy_K[i,j] * tab_S2_K[i,j]
            tab_conv_Byy_S2_K[i,j] = TAB_Byy_K[i,j] * tab_S2_K[i,j]
            #
            tab_conv_Bxx_dS2dx_K[i,j] = TAB_Bxx_K[i,j] * tab_dS2dx_K[i,j]
            tab_conv_Bxy_dS2dx_K[i,j] = TAB_Bxy_K[i,j] * tab_dS2dx_K[i,j]
            tab_conv_Bxy_dS2dy_K[i,j] = TAB_Bxy_K[i,j] * tab_dS2dy_K[i,j]
            tab_conv_Byy_dS2dy_K[i,j] = TAB_Byy_K[i,j] * tab_dS2dy_K[i,j]
            ####
        end
    end    
end

##################################################
# Compute the inverse FFT of convolution results in L_v S2
##################################################
function compute_inv_FFT_conv_S2!(cnt_FFT_S2)
    tab_conv_Bxx_S2_K = cnt_FFT_S2.tab_conv_Bxx_S2_K
    tab_conv_Bxy_S2_K = cnt_FFT_S2.tab_conv_Bxy_S2_K
    tab_conv_Byy_S2_K = cnt_FFT_S2.tab_conv_Byy_S2_K
    #
    tab_conv_Bxx_S2_R = cnt_FFT_S2.tab_conv_Bxx_S2_R
    tab_conv_Bxy_S2_R = cnt_FFT_S2.tab_conv_Bxy_S2_R
    tab_conv_Byy_S2_R = cnt_FFT_S2.tab_conv_Byy_S2_R
    #
    tab_conv_Bxx_dS2dx_K = cnt_FFT_S2.tab_conv_Bxx_dS2dx_K
    tab_conv_Bxy_dS2dx_K = cnt_FFT_S2.tab_conv_Bxy_dS2dx_K
    tab_conv_Bxy_dS2dy_K = cnt_FFT_S2.tab_conv_Bxy_dS2dy_K
    tab_conv_Byy_dS2dy_K = cnt_FFT_S2.tab_conv_Byy_dS2dy_K
    #
    tab_conv_Bxx_dS2dx_R = cnt_FFT_S2.tab_conv_Bxx_dS2dx_R
    tab_conv_Bxy_dS2dx_R = cnt_FFT_S2.tab_conv_Bxy_dS2dx_R
    tab_conv_Bxy_dS2dy_R = cnt_FFT_S2.tab_conv_Bxy_dS2dy_R
    tab_conv_Byy_dS2dy_R = cnt_FFT_S2.tab_conv_Byy_dS2dy_R 

    # Compute the inverse FFT of the convolution
    mul!(tab_conv_Bxx_S2_R, PLAN_K_TO_R, tab_conv_Bxx_S2_K)
    mul!(tab_conv_Bxy_S2_R, PLAN_K_TO_R, tab_conv_Bxy_S2_K)
    mul!(tab_conv_Byy_S2_R, PLAN_K_TO_R, tab_conv_Byy_S2_K)
    #
    mul!(tab_conv_Bxx_dS2dx_R, PLAN_K_TO_R, tab_conv_Bxx_dS2dx_K)
    mul!(tab_conv_Bxy_dS2dx_R, PLAN_K_TO_R, tab_conv_Bxy_dS2dx_K)
    mul!(tab_conv_Bxy_dS2dy_R, PLAN_K_TO_R, tab_conv_Bxy_dS2dy_K)
    mul!(tab_conv_Byy_dS2dy_R, PLAN_K_TO_R, tab_conv_Byy_dS2dy_K)
end


##################################################
# Compute_Lv using FFTs
# + Threads.@threads :static for l2
# + WITHOUT allocation
##################################################
function compute_Lv_contribution!(TAB_F, TAB_dF_dx, TAB_dF_dy, indx)
    
    # Loop over v2 index 
    Threads.@threads :static for l2 in 1:NV2
        #####
        tid = Threads.threadid()
        #####
        # Allocate thread-local temporary arrays to avoid data races
        tab_S2_slice     = TAB_S2_SLICE_PARALLEL[tid]
        tab_dS2_dx_slice = TAB_DS2_DX_PARALLEL[tid]
        tab_dS2_dy_slice = TAB_DS2_DY_PARALLEL[tid]
        Integrand_x      = TAB_Ix_PARALLEL[tid]
        Integrand_y      = TAB_Iy_PARALLEL[tid]
        div_Int          = TAB_DIV_I_PARALLEL[tid]

        #####
        # From here, v2 = l2 is fixed.
        # Build S2_slice(v1) = S2(v1, v2) if indx == 0 or
        #                    = S2(v2, v1) if indx == 1.
        slice_2D_field(indx, l2, tab_S2_slice)

        #####
        # Compute gradients with respect to v1 for S2_slice
        compute_1D_field_gradients!(tab_S2_slice, tab_dS2_dx_slice, tab_dS2_dy_slice)
        
        #####
        # FFT part
        #####
        cnt_S2_loc = CONTAINER_FFT_S2[tid]

        # Prepare arrays in the R-space for this slice of S2
        prepare_fields_S2!(tab_S2_slice, tab_dS2_dx_slice, tab_dS2_dy_slice, cnt_S2_loc)
        # Compute FFT of S2 and p_v S2
        compute_FFT_S2!(cnt_S2_loc)
        # Compute all the convolutions
        compute_convolution_S2!(cnt_S2_loc)
        # Compute the inv FFT of all the convolutions
        compute_inv_FFT_conv_S2!(cnt_S2_loc) # This functions is the most expensive

        # Fetch pointers locally
        tab_conv_Bxx_S2 = cnt_S2_loc.tab_conv_Bxx_S2_R
        tab_conv_Bxy_S2 = cnt_S2_loc.tab_conv_Bxy_S2_R
        tab_conv_Byy_S2 = cnt_S2_loc.tab_conv_Byy_S2_R
        #
        tab_conv_Bxx_dS2dx = cnt_S2_loc.tab_conv_Bxx_dS2dx_R
        tab_conv_Bxy_dS2dx = cnt_S2_loc.tab_conv_Bxy_dS2dx_R
        tab_conv_Bxy_dS2dy = cnt_S2_loc.tab_conv_Bxy_dS2dy_R
        tab_conv_Byy_dS2dy = cnt_S2_loc.tab_conv_Byy_dS2dy_R

        #####
        # For each v1 = l1, assemble the integrand
        #####
        @inbounds @simd for l1 in 1:NV2
            i, j = TAB_i[l1], TAB_j[l1]

            # F and p_v F
            Fv    = TAB_F[l1]
            dF_dx = TAB_dF_dx[l1]
            dF_dy = TAB_dF_dy[l1]
            
            # S2 and p_v S2
            S2_v  = tab_S2_slice[l1] 
            dS2_x = tab_dS2_dx_slice[l1]
            dS2_y = tab_dS2_dy_slice[l1]

            # Convolutions terms in the real space
            int_Bxx_S2   = tab_conv_Bxx_S2[i,j]
            int_Bxy_S2   = tab_conv_Bxy_S2[i,j]
            int_Byy_S2   = tab_conv_Byy_S2[i,j]
            #
            int_Bxx_dS2x = tab_conv_Bxx_dS2dx[i,j]
            int_Bxy_dS2x = tab_conv_Bxy_dS2dx[i,j]
            int_Bxy_dS2y = tab_conv_Bxy_dS2dy[i,j]
            int_Byy_dS2y = tab_conv_Byy_dS2dy[i,j]

            #####
            # Groups all the terms
            #####
            term_xx = (dF_dx * int_Bxx_S2) - (S2_v * TAB_INT_Bxx_dFx[l1]) + (dS2_x * TAB_INT_Bxx_F[l1]) - (Fv * int_Bxx_dS2x)
            term_xy = (dF_dy * int_Bxy_S2) - (S2_v * TAB_INT_Bxy_dFy[l1]) + (dS2_y * TAB_INT_Bxy_F[l1]) - (Fv * int_Bxy_dS2y)
            #
            term_yx = (dF_dx * int_Bxy_S2) - (S2_v * TAB_INT_Bxy_dFx[l1]) + (dS2_x * TAB_INT_Bxy_F[l1]) - (Fv * int_Bxy_dS2x)
            term_yy = (dF_dy * int_Byy_S2) - (S2_v * TAB_INT_Byy_dFy[l1]) + (dS2_y * TAB_INT_Byy_F[l1]) - (Fv * int_Byy_dS2y)

            # Accumulators for the integral result
            Integrand_x[l1] = (term_xx + term_xy) * DV2
            Integrand_y[l1] = (term_yx + term_yy) * DV2
        end

        #####
        # Compute the divergence of TAB_Integrand (p_v1x Integrand_x + p_v1y Integrand_y)
        # by finite difference.
        # Gives p_tS contribution for this v2 = l2
        compute_1D_field_div!(Integrand_x, Integrand_y, div_Int)
        
        #####
        # Accumulate this v2 contribution into the global array.
        # ATTENTION, the accumulation depends on the direction = v1 or v2
        if indx == 0     # We are computing L_v1 S_2(v1, v2)
            @simd for l1 in 1:NV2
                TAB_dS2_dt[l1, l2] += div_Int[l1]
            end
        #
        elseif indx == 1 # We are computing L_v2 S_2(v1, v2)
            @simd for l1 in 1:NV2
                TAB_dS2_dt[l2, l1] += div_Int[l1]
            end
        end
    end
end
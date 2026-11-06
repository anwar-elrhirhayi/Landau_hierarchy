####################################################################################################
# PDE for S(v1, v2,t) for Landau - finite differences
####################################################################################################

####################################################################################################
# For a given TAB_S2, compute contribution L_v1 S(v1, v2, v3) if indx == 0
#                                          L_v2 S(v1, v2, v3) if indx == 1
#                                          L_v3 S(v1, v2, v3) if indx == 2
####################################################################################################

##################################################
# Containers for the parallelisation in L_v S3
##################################################
struct container_FFT_S3
    # Arrays in the R-space
    tab_S3_R::Array{Float64, 2}
    tab_dS3dx_R::Array{Float64, 2}
    tab_dS3dy_R::Array{Float64, 2}
    #
    tab_conv_Bxx_S3_R::Array{Float64, 2}
    tab_conv_Bxy_S3_R::Array{Float64, 2}
    tab_conv_Byy_S3_R::Array{Float64, 2}
    #
    tab_conv_Bxx_dS3dx_R::Array{Float64, 2}
    tab_conv_Bxy_dS3dx_R::Array{Float64, 2}
    tab_conv_Bxy_dS3dy_R::Array{Float64, 2}
    tab_conv_Byy_dS3dy_R::Array{Float64, 2}
    
    # Arrays in the K-space
    tab_S3_K::Array{Complex{Float64}, 2}
    tab_dS3dx_K::Array{Complex{Float64}, 2}
    tab_dS3dy_K::Array{Complex{Float64}, 2}
    #
    tab_conv_Bxx_S3_K::Array{Complex{Float64}, 2}
    tab_conv_Bxy_S3_K::Array{Complex{Float64}, 2}
    tab_conv_Byy_S3_K::Array{Complex{Float64}, 2}
    #
    tab_conv_Bxx_dS3dx_K::Array{Complex{Float64}, 2}
    tab_conv_Bxy_dS3dx_K::Array{Complex{Float64}, 2}
    tab_conv_Bxy_dS3dy_K::Array{Complex{Float64}, 2}
    tab_conv_Byy_dS3dy_K::Array{Complex{Float64}, 2}
end

const CONTAINER_FFT_S3 = [ container_FFT_S3(zeros(Float64,NV_DOUBLE,NV_DOUBLE),
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
# Prepare the FFT of S3, p_v S3 
##################################################
function prepare_fields_S3!(tab_S3_slice, tab_dS3_dx_slice, tab_dS3_dy_slice, cnt_FFT_S3)
    #
    tab_S3_R     = cnt_FFT_S3.tab_S3_R
    tab_dS3dx_R  = cnt_FFT_S3.tab_dS3dx_R
    tab_dS3dy_R  = cnt_FFT_S3.tab_dS3dy_R
    #
    fill!(tab_S3_R   , 0.0)
    fill!(tab_dS3dx_R, 0.0)
    fill!(tab_dS3dy_R, 0.0)
    #
    @inbounds for l = 1:NV2
        #####
        i, j = TAB_i[l], TAB_j[l]
        #
        tab_S3_R[i,j]    = tab_S3_slice[l]
        tab_dS3dx_R[i,j] = tab_dS3_dx_slice[l]
        tab_dS3dy_R[i,j] = tab_dS3_dy_slice[l]
        #####
    end
end

##################################################
# Compute the FFT of S3 slices and its gradients
##################################################
function compute_FFT_S3!(cnt_FFT_S3)
    tab_S3_R        = cnt_FFT_S3.tab_S3_R
    tab_dS3dx_R     = cnt_FFT_S3.tab_dS3dx_R
    tab_dS3dy_R     = cnt_FFT_S3.tab_dS3dy_R
    #
    tab_S3_K        = cnt_FFT_S3.tab_S3_K
    tab_dS3dx_K     = cnt_FFT_S3.tab_dS3dx_K
    tab_dS3dy_K     = cnt_FFT_S3.tab_dS3dy_K

    # Compute FFT
    mul!(tab_S3_K,    PLAN_R_TO_K, tab_S3_R)
    mul!(tab_dS3dx_K, PLAN_R_TO_K, tab_dS3dx_R)
    mul!(tab_dS3dy_K, PLAN_R_TO_K, tab_dS3dy_R)
end

##################################################
# Compute the convolution in L_v S3
##################################################
function compute_convolution_S3!(cnt_FFT_S3)
    tab_S3_K        = cnt_FFT_S3.tab_S3_K
    tab_dS3dx_K     = cnt_FFT_S3.tab_dS3dx_K
    tab_dS3dy_K     = cnt_FFT_S3.tab_dS3dy_K
    #
    tab_conv_Bxx_S3_K = cnt_FFT_S3.tab_conv_Bxx_S3_K
    tab_conv_Bxy_S3_K = cnt_FFT_S3.tab_conv_Bxy_S3_K
    tab_conv_Byy_S3_K = cnt_FFT_S3.tab_conv_Byy_S3_K
    #
    tab_conv_Bxx_dS3dx_K = cnt_FFT_S3.tab_conv_Bxx_dS3dx_K
    tab_conv_Bxy_dS3dx_K = cnt_FFT_S3.tab_conv_Bxy_dS3dx_K
    tab_conv_Bxy_dS3dy_K = cnt_FFT_S3.tab_conv_Bxy_dS3dy_K
    tab_conv_Byy_dS3dy_K = cnt_FFT_S3.tab_conv_Byy_dS3dy_K

    # Compute the convolution
    @inbounds for j = 1:NV_DOUBLE
        @simd for i = 1:NV+1
            ####
            s3_k   = tab_S3_K[i,j]
            ds3x_k = tab_dS3dx_K[i,j]
            ds3y_k = tab_dS3dy_K[i,j]
            #
            bxx = TAB_Bxx_K[i,j]
            bxy = TAB_Bxy_K[i,j]
            byy = TAB_Byy_K[i,j]
            ####
            tab_conv_Bxx_S3_K[i,j] = bxx * s3_k
            tab_conv_Bxy_S3_K[i,j] = bxy * s3_k
            tab_conv_Byy_S3_K[i,j] = byy * s3_k
            #
            tab_conv_Bxx_dS3dx_K[i,j] = bxx * ds3x_k
            tab_conv_Bxy_dS3dx_K[i,j] = bxy * ds3x_k
            tab_conv_Bxy_dS3dy_K[i,j] = bxy * ds3y_k
            tab_conv_Byy_dS3dy_K[i,j] = byy * ds3y_k
            ####
        end
    end    
end

##################################################
# Compute the inverse FFT of convolution results in L_v S3
##################################################
function compute_inv_FFT_conv_S3!(cnt_FFT_S3)
    tab_conv_Bxx_S3_K = cnt_FFT_S3.tab_conv_Bxx_S3_K
    tab_conv_Bxy_S3_K = cnt_FFT_S3.tab_conv_Bxy_S3_K
    tab_conv_Byy_S3_K = cnt_FFT_S3.tab_conv_Byy_S3_K
    #
    tab_conv_Bxx_S3_R = cnt_FFT_S3.tab_conv_Bxx_S3_R
    tab_conv_Bxy_S3_R = cnt_FFT_S3.tab_conv_Bxy_S3_R
    tab_conv_Byy_S3_R = cnt_FFT_S3.tab_conv_Byy_S3_R
    #
    tab_conv_Bxx_dS3dx_K = cnt_FFT_S3.tab_conv_Bxx_dS3dx_K
    tab_conv_Bxy_dS3dx_K = cnt_FFT_S3.tab_conv_Bxy_dS3dx_K
    tab_conv_Bxy_dS3dy_K = cnt_FFT_S3.tab_conv_Bxy_dS3dy_K
    tab_conv_Byy_dS3dy_K = cnt_FFT_S3.tab_conv_Byy_dS3dy_K
    #
    tab_conv_Bxx_dS3dx_R = cnt_FFT_S3.tab_conv_Bxx_dS3dx_R
    tab_conv_Bxy_dS3dx_R = cnt_FFT_S3.tab_conv_Bxy_dS3dx_R
    tab_conv_Bxy_dS3dy_R = cnt_FFT_S3.tab_conv_Bxy_dS3dy_R
    tab_conv_Byy_dS3dy_R = cnt_FFT_S3.tab_conv_Byy_dS3dy_R 

    # Compute the inverse FFT of the convolution
    mul!(tab_conv_Bxx_S3_R, PLAN_K_TO_R, tab_conv_Bxx_S3_K)
    mul!(tab_conv_Bxy_S3_R, PLAN_K_TO_R, tab_conv_Bxy_S3_K)
    mul!(tab_conv_Byy_S3_R, PLAN_K_TO_R, tab_conv_Byy_S3_K)
    #
    mul!(tab_conv_Bxx_dS3dx_R, PLAN_K_TO_R, tab_conv_Bxx_dS3dx_K)
    mul!(tab_conv_Bxy_dS3dx_R, PLAN_K_TO_R, tab_conv_Bxy_dS3dx_K)
    mul!(tab_conv_Bxy_dS3dy_R, PLAN_K_TO_R, tab_conv_Bxy_dS3dy_K)
    mul!(tab_conv_Byy_dS3dy_R, PLAN_K_TO_R, tab_conv_Byy_dS3dy_K)
end

##################################################
# Compute Lv
##################################################
function compute_Lv_contribution_S3!(TAB_F, TAB_dF_dx, TAB_dF_dy, indx)
    # Loop over v2 and v3 --> index l2, l3

    Threads.@threads :static for l2_l3 in CartesianIndices((1:NV2, 1:NV2))
        l2, l3 = l2_l3[1], l2_l3[2]
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
        # FFT part
        #####
        cnt_S3_loc = CONTAINER_FFT_S3[tid]

        # Prepare arrays in the R-space for this slice of S3
        prepare_fields_S3!(tab_S3_slice, tab_dS3_dx_slice, tab_dS3_dy_slice, cnt_S3_loc)
        # Compute FFT of S3 and p_v S3
        compute_FFT_S3!(cnt_S3_loc)
        # Compute all the convolutions
        compute_convolution_S3!(cnt_S3_loc)
        # Compute the inv FFT of all the convolutions
        compute_inv_FFT_conv_S3!(cnt_S3_loc) # This functions is the most expensive

        tab_conv_Bxx_S3 = cnt_S3_loc.tab_conv_Bxx_S3_R
        tab_conv_Bxy_S3 = cnt_S3_loc.tab_conv_Bxy_S3_R
        tab_conv_Byy_S3 = cnt_S3_loc.tab_conv_Byy_S3_R
        #
        tab_conv_Bxx_dS3dx = cnt_S3_loc.tab_conv_Bxx_dS3dx_R
        tab_conv_Bxy_dS3dx = cnt_S3_loc.tab_conv_Bxy_dS3dx_R
        tab_conv_Bxy_dS3dy = cnt_S3_loc.tab_conv_Bxy_dS3dy_R
        tab_conv_Byy_dS3dy = cnt_S3_loc.tab_conv_Byy_dS3dy_R

        # For each v1 = l1, we compute the integral over v = l
        @inbounds @simd for l1=1:NV2
            i, j = TAB_i[l1], TAB_j[l1]

            # F and p_v F
            Fv    = TAB_F[l1]
            dF_dx = TAB_dF_dx[l1]
            dF_dy = TAB_dF_dy[l1]
            
            # S3 and p_v S3
            S3_v  = tab_S3_slice[l1] 
            dS3_x = tab_dS3_dx_slice[l1]
            dS3_y = tab_dS3_dy_slice[l1]

            # Convolutions terms in the real space
            int_Bxx_S3   = tab_conv_Bxx_S3[i,j]
            int_Bxy_S3   = tab_conv_Bxy_S3[i,j]
            int_Byy_S3   = tab_conv_Byy_S3[i,j]
            #
            int_Bxx_dS3x = tab_conv_Bxx_dS3dx[i,j]
            int_Bxy_dS3x = tab_conv_Bxy_dS3dx[i,j]
            int_Bxy_dS3y = tab_conv_Bxy_dS3dy[i,j]
            int_Byy_dS3y = tab_conv_Byy_dS3dy[i,j]

            # Compute the components of the integral in the Flux term
            term_xx = (dF_dx * int_Bxx_S3 * DV2) - S3_v * TAB_INT_Bxx_dFx[l1] + dS3_x * TAB_INT_Bxx_F[l1] - (Fv * int_Bxx_dS3x * DV2)
            term_xy = (dF_dy * int_Bxy_S3 * DV2) - S3_v * TAB_INT_Bxy_dFy[l1] + dS3_y * TAB_INT_Bxy_F[l1] - (Fv * int_Bxy_dS3y * DV2)
            # 
            term_yx = (dF_dx * int_Bxy_S3 * DV2) - S3_v * TAB_INT_Bxy_dFx[l1] + dS3_x * TAB_INT_Bxy_F[l1] - (Fv * int_Bxy_dS3x * DV2)
            term_yy = (dF_dy * int_Byy_S3 * DV2) - S3_v * TAB_INT_Byy_dFy[l1] + dS3_y * TAB_INT_Byy_F[l1] - (Fv * int_Byy_dS3y * DV2)
            
            # Accumulators for the integral result
            # ATTENTION to the factor 1/2
            Integrand_x[l1] = 0.5 * (term_xx + term_xy)
            Integrand_y[l1] = 0.5 * (term_yx + term_yy)
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
            @simd for l1 in 1:NV2
                TAB_dS3_dt[l1, l2, l3] += div_Int[l1]
            end
        #
        elseif indx == 1 # We are computing L_v2 S_3(v3, v1, v2)
            @simd for l1 in 1:NV2
                TAB_dS3_dt[l3, l1, l2] += div_Int[l1]
            end
        #    
        elseif indx == 2 # We are computing L_v3 S_3(v2, v3, v1)
            @simd for l1 in 1:NV2
                TAB_dS3_dt[l2, l3, l1] += div_Int[l1]
            end
        end
    end
end


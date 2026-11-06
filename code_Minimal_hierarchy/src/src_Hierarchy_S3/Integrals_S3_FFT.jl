##################################################
# Arrays for the S2 convolutions used by S3
##################################################
const TAB_INT_I2_xx = zeros(Float64, NV2, NV2)
const TAB_INT_I2_xy = zeros(Float64, NV2, NV2)
const TAB_INT_I2_yy = zeros(Float64, NV2, NV2)

const TAB_INT_I3_x = zeros(Float64, NV2, NV2)
const TAB_INT_I3_y = zeros(Float64, NV2, NV2)

const TAB_INT_S2_Bxx = zeros(Float64, NV2, NV2)
const TAB_INT_S2_Bxy = zeros(Float64, NV2, NV2)
const TAB_INT_S2_Byy = zeros(Float64, NV2, NV2)

const TAB_INTEGRALS_CONV_R = zeros(Float64, NV_DOUBLE, NV_DOUBLE)

##################################################
# Compute convolution
##################################################
function compute_integral_convolution!(conv_R, conv_K, kernel_K, field_K)
    @inbounds @simd for index in eachindex(conv_K, kernel_K, field_K)
        conv_K[index] = kernel_K[index] * field_K[index]
    end
    mul!(conv_R, PLAN_K_TO_R, conv_K)
end

##################################################
# Copy the convolution back to the physical velocity grid
##################################################
function copy_convolution_to_grid!(destination, convolution)
    @inbounds @simd for l in 1:NV2
        destination[l] = convolution[TAB_i[l], TAB_j[l]] * DV2
    end
end

####################################################################################################
# Precompute the F integrals used by the Lv operators for S2 and S3
####################################################################################################
function precompute_F_B_integrals!()
    prepare_arrays_F!()
    FFT_fields_F!()

    compute_integral_convolution!(TAB_INTEGRALS_CONV_R, CONV_K, TAB_Bxx_K, TAB_F_K)
    copy_convolution_to_grid!(TAB_INT_Bxx_F, TAB_INTEGRALS_CONV_R)

    compute_integral_convolution!(TAB_INTEGRALS_CONV_R, CONV_K, TAB_Bxy_K, TAB_F_K)
    copy_convolution_to_grid!(TAB_INT_Bxy_F, TAB_INTEGRALS_CONV_R)

    compute_integral_convolution!(TAB_INTEGRALS_CONV_R, CONV_K, TAB_Byy_K, TAB_F_K)
    copy_convolution_to_grid!(TAB_INT_Byy_F, TAB_INTEGRALS_CONV_R)

    compute_integral_convolution!(TAB_INTEGRALS_CONV_R, CONV_K, TAB_Bxx_K, TAB_dFdx_K)
    copy_convolution_to_grid!(TAB_INT_Bxx_dFx, TAB_INTEGRALS_CONV_R)

    compute_integral_convolution!(TAB_INTEGRALS_CONV_R, CONV_K, TAB_Bxy_K, TAB_dFdx_K)
    copy_convolution_to_grid!(TAB_INT_Bxy_dFx, TAB_INTEGRALS_CONV_R)

    compute_integral_convolution!(TAB_INTEGRALS_CONV_R, CONV_K, TAB_Bxy_K, TAB_dFdy_K)
    copy_convolution_to_grid!(TAB_INT_Bxy_dFy, TAB_INTEGRALS_CONV_R)

    compute_integral_convolution!(TAB_INTEGRALS_CONV_R, CONV_K, TAB_Byy_K, TAB_dFdy_K)
    copy_convolution_to_grid!(TAB_INT_Byy_dFy, TAB_INTEGRALS_CONV_R)
end

##################################################
# Compute (v1, v2) \mapsto int dv S2(v1, v) B(v, v2).
##################################################
function precompute_S2_B_integral!()
    fill!(TAB_INT_S2_Bxx, 0.0)
    fill!(TAB_INT_S2_Bxy, 0.0)
    fill!(TAB_INT_S2_Byy, 0.0)

    Threads.@threads :static for l1 in 1:NV2
        tid = Threads.threadid()
        container = CONTAINER_FFT_S3[tid]

        S2_R_loc = container.tab_S3_R
        S2_K_loc = container.tab_S3_K
        conv_R_loc = container.tab_conv_Bxx_S3_R
        conv_K_loc = container.tab_conv_Bxx_S3_K

        fill!(S2_R_loc, 0.0)
        @inbounds for l in 1:NV2
            S2_R_loc[TAB_i[l], TAB_j[l]] = TAB_S2[l1, l]
        end
        mul!(S2_K_loc, PLAN_R_TO_K, S2_R_loc)

        compute_integral_convolution!(conv_R_loc, conv_K_loc, TAB_Bxx_K, S2_K_loc)
        @inbounds @simd for l2 in 1:NV2
            TAB_INT_S2_Bxx[l1, l2] = conv_R_loc[TAB_i[l2], TAB_j[l2]] * DV2
        end

        compute_integral_convolution!(conv_R_loc, conv_K_loc, TAB_Bxy_K, S2_K_loc)
        @inbounds @simd for l2 in 1:NV2
            TAB_INT_S2_Bxy[l1, l2] = conv_R_loc[TAB_i[l2], TAB_j[l2]] * DV2
        end

        compute_integral_convolution!(conv_R_loc, conv_K_loc, TAB_Byy_K, S2_K_loc)
        @inbounds @simd for l2 in 1:NV2
            TAB_INT_S2_Byy[l1, l2] = conv_R_loc[TAB_i[l2], TAB_j[l2]] * DV2
        end
    end
end

##################################################
# Compute the I2 and I3 terms used by the NL S3 term.
##################################################
function TAB_INTEGRALS_NL!()
    fill!(TAB_INT_I2_xx, 0.0)
    fill!(TAB_INT_I2_xy, 0.0)
    fill!(TAB_INT_I2_yy, 0.0)
    fill!(TAB_INT_I3_x, 0.0)
    fill!(TAB_INT_I3_y, 0.0)

    Threads.@threads :static for l2 in 1:NV2
        tid = Threads.threadid()
        container = CONTAINER_FFT_S3[tid]

        S2_R_loc = container.tab_S3_R
        dS2dx_R_loc = container.tab_dS3dx_R
        dS2dy_R_loc = container.tab_dS3dy_R
        S2_K_loc = container.tab_S3_K
        dS2dx_K_loc = container.tab_dS3dx_K
        dS2dy_K_loc = container.tab_dS3dy_K
        conv_R_loc = container.tab_conv_Bxx_S3_R
        conv_K_loc = container.tab_conv_Bxx_S3_K

        fill!(S2_R_loc, 0.0)
        fill!(dS2dx_R_loc, 0.0)
        fill!(dS2dy_R_loc, 0.0)

        @inbounds for l in 1:NV2
            i, j = TAB_i[l], TAB_j[l]
            S2_R_loc[i, j] = TAB_S2[l, l2]
            dS2dx_R_loc[i, j] = TAB_dS2_dx[l, l2]
            dS2dy_R_loc[i, j] = TAB_dS2_dy[l, l2]
        end

        mul!(S2_K_loc, PLAN_R_TO_K, S2_R_loc)
        mul!(dS2dx_K_loc, PLAN_R_TO_K, dS2dx_R_loc)
        mul!(dS2dy_K_loc, PLAN_R_TO_K, dS2dy_R_loc)

        # I2 = B * S2
        compute_integral_convolution!(conv_R_loc, conv_K_loc, TAB_Bxx_K, S2_K_loc)
        @inbounds @simd for l1 in 1:NV2
            TAB_INT_I2_xx[l1, l2] = conv_R_loc[TAB_i[l1], TAB_j[l1]] * DV2
        end

        compute_integral_convolution!(conv_R_loc, conv_K_loc, TAB_Bxy_K, S2_K_loc)
        @inbounds @simd for l1 in 1:NV2
            TAB_INT_I2_xy[l1, l2] = conv_R_loc[TAB_i[l1], TAB_j[l1]] * DV2
        end

        compute_integral_convolution!(conv_R_loc, conv_K_loc, TAB_Byy_K, S2_K_loc)
        @inbounds @simd for l1 in 1:NV2
            TAB_INT_I2_yy[l1, l2] = conv_R_loc[TAB_i[l1], TAB_j[l1]] * DV2
        end

        # I3_x = Bxx * dS2_dx + Bxy * dS2_dy
        compute_integral_convolution!(conv_R_loc, conv_K_loc, TAB_Bxx_K, dS2dx_K_loc)
        @inbounds @simd for l1 in 1:NV2
            TAB_INT_I3_x[l1, l2] = conv_R_loc[TAB_i[l1], TAB_j[l1]] * DV2
        end

        compute_integral_convolution!(conv_R_loc, conv_K_loc, TAB_Bxy_K, dS2dy_K_loc)
        @inbounds @simd for l1 in 1:NV2
            TAB_INT_I3_x[l1, l2] += conv_R_loc[TAB_i[l1], TAB_j[l1]] * DV2
        end

        # I3_y = Bxy * dS2_dx + Byy * dS2_dy
        compute_integral_convolution!(conv_R_loc, conv_K_loc, TAB_Bxy_K, dS2dx_K_loc)
        @inbounds @simd for l1 in 1:NV2
            TAB_INT_I3_y[l1, l2] = conv_R_loc[TAB_i[l1], TAB_j[l1]] * DV2
        end

        compute_integral_convolution!(conv_R_loc, conv_K_loc, TAB_Byy_K, dS2dy_K_loc)
        @inbounds @simd for l1 in 1:NV2
            TAB_INT_I3_y[l1, l2] += conv_R_loc[TAB_i[l1], TAB_j[l1]] * DV2
        end
    end
end

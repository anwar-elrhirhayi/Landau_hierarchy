##################################################
# Landau equation in 2D
# Discretization of velocity space with double integral
# Explicit Euler scheme for time integration
##################################################
# Compute the flux, Flux(v):
# p_t F = div Flux(v)
# Flux(v) = p_v int dv' B(v, v') [- F(v) p_v'F + F(v') p_vF]
# Implementation: via FFT
# Results stored in TAB_Fluxx, TAB_Fluxy
##################################################


##################################################
# Compute the total flux, by computing each convolution
# one-by-one
##################################################
function TAB_Flux!()

    # Prepare the FT of F and p_v F
    prepare_arrays_F!()
    FFT_fields_F!()
    
    ########
    # Compute each convolution one-by-one
    ########

    # Term Bxx * F
    compute_convolution!(TAB_Bxx_K, TAB_F_K)
    compute_inv_FFT!(CONV_R_Bxx_F, CONV_K)

    # Term Bxx * dFdx
    compute_convolution!(TAB_Bxx_K, TAB_dFdx_K)
    compute_inv_FFT!(CONV_R_Bxx_dFdx, CONV_K)

    # Term Bxy * F
    compute_convolution!(TAB_Bxy_K, TAB_F_K)
    compute_inv_FFT!(CONV_R_Bxy_F, CONV_K)

    # Term Bxy * dFdy
    compute_convolution!(TAB_Bxy_K, TAB_dFdy_K)
    compute_inv_FFT!(CONV_R_Bxy_dFdy, CONV_K)

    # Term Bxy * dFdx
    compute_convolution!(TAB_Bxy_K, TAB_dFdx_K) 
    compute_inv_FFT!(CONV_R_Bxy_dFdx, CONV_K)
    
    # Term Byy * F
    compute_convolution!(TAB_Byy_K, TAB_F_K)
    compute_inv_FFT!(CONV_R_Byy_F, CONV_K)

    # Term Byy * dFdy
    compute_convolution!(TAB_Byy_K, TAB_dFdy_K)
    compute_inv_FFT!(CONV_R_Byy_dFdy, CONV_K)
    
    @inbounds @simd for l=1:NV2
        #####
        i, j = TAB_i[l], TAB_j[l]
        
        # Read F, p_v F
        Fv     = TAB_F[l]
        dFdx_v = TAB_dF_dx[l]
        dFdy_v = TAB_dF_dy[l]

        # Compute flux along the x-direction
        term_xx = dFdx_v * CONV_R_Bxx_F[i,j] - Fv * CONV_R_Bxx_dFdx[i,j]
        term_xy = dFdy_v * CONV_R_Bxy_F[i,j] - Fv * CONV_R_Bxy_dFdy[i,j]
        TAB_Fluxx[l] = (term_xx + term_xy) * 0.5 * DV^2

        # Compute flux along the y-direction
        term_yx = dFdx_v * CONV_R_Bxy_F[i,j] - Fv * CONV_R_Bxy_dFdx[i,j]
        term_yy = dFdy_v * CONV_R_Byy_F[i,j] - Fv * CONV_R_Byy_dFdy[i,j]
        TAB_Fluxy[l] = (term_yx + term_yy) * 0.5 * DV^2
        #####
    end
end

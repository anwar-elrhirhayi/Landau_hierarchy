##################################################
# General FFT arrays
##################################################
const NV_DOUBLE = NV * 2  # Used for the FFT
#
const TAB_F_R    = zeros(Float64, NV_DOUBLE, NV_DOUBLE)
const TAB_F_K    = zeros(Complex{Float64}, NV+1, NV_DOUBLE)
#
const PLAN_FFT_R_TO_K = plan_rfft(TAB_F_R, flags=FFTW.MEASURE)
const PLAN_FFT_K_TO_R = plan_irfft(TAB_F_K, NV_DOUBLE, flags=FFTW.MEASURE)


const PLAN_R_TO_K =  plan_rfft(TAB_F_R           , flags=FFTW.MEASURE)
const PLAN_K_TO_R = plan_irfft(TAB_F_K, NV_DOUBLE, flags=FFTW.MEASURE) 
##################################################
# Arrays used for the FFT of p_v F
##################################################
const TAB_dFdx_R = zeros(Float64, NV_DOUBLE, NV_DOUBLE)
const TAB_dFdy_R = zeros(Float64, NV_DOUBLE, NV_DOUBLE)
#
const TAB_dFdx_K = zeros(Complex{Float64}, NV+1, NV_DOUBLE)
const TAB_dFdy_K = zeros(Complex{Float64}, NV+1, NV_DOUBLE)

##################################################
# Arrays used for the convolution
##################################################
const CONV_K = zeros(Complex{Float64}, NV+1, NV_DOUBLE)
const CONV_R_Bxx_F    = zeros(Float64, NV_DOUBLE, NV_DOUBLE)
const CONV_R_Bxx_dFdx = zeros(Float64, NV_DOUBLE, NV_DOUBLE) 
const CONV_R_Bxy_F    = zeros(Float64, NV_DOUBLE, NV_DOUBLE)
const CONV_R_Bxy_dFdy = zeros(Float64, NV_DOUBLE, NV_DOUBLE)
const CONV_R_Bxy_dFdx = zeros(Float64, NV_DOUBLE, NV_DOUBLE)
const CONV_R_Byy_F    = zeros(Float64, NV_DOUBLE, NV_DOUBLE)
const CONV_R_Byy_dFdy = zeros(Float64, NV_DOUBLE, NV_DOUBLE)

##################################################
# Prepare the arrays of F, p_v F for the FFT
##################################################
function prepare_arrays_F!()
    fill!(TAB_F_R, 0.0)
    fill!(TAB_dFdx_R, 0.0)
    fill!(TAB_dFdy_R, 0.0)
    # 
    @inbounds for l = 1:NV2
        i, j = TAB_i[l], TAB_j[l]
        TAB_F_R[i,j]    = TAB_F[l]
        TAB_dFdx_R[i,j] = TAB_dF_dx[l]
        TAB_dFdy_R[i,j] = TAB_dF_dy[l]
    end
end

##################################################
# Compute FFT of F and p_v F
##################################################
function FFT_fields_F!()
    mul!(TAB_F_K,    PLAN_FFT_R_TO_K, TAB_F_R)
    mul!(TAB_dFdx_K, PLAN_FFT_R_TO_K, TAB_dFdx_R)
    mul!(TAB_dFdy_K, PLAN_FFT_R_TO_K, TAB_dFdy_R)
end

##################################################
# Compute the convolution between Kernel_K and Field_K
# The result is stored in CONV_K
##################################################
function compute_convolution!(Kernel_K, Field_K)
    @inbounds for j = 1:NV_DOUBLE
        @simd for i = 1:NV+1
            CONV_K[i,j] = Kernel_K[i,j] * Field_K[i,j]
        end
    end
end

##################################################
# Compute the inverse FFT of a field Tab_k
################################################## 
function compute_inv_FFT!(Tab_R, Tab_K)
    # Inverse FFT
    mul!(Tab_R, PLAN_FFT_K_TO_R, Tab_K)
end

##################################################
# Arrays for B
##################################################
const TAB_Bxx_K = zeros(Complex{Float64}, NV+1, NV_DOUBLE)
const TAB_Bxy_K = zeros(Complex{Float64}, NV+1, NV_DOUBLE)
const TAB_Byy_K = zeros(Complex{Float64}, NV+1, NV_DOUBLE)
#
const TAB_Bxx_R = zeros(Float64, NV_DOUBLE, NV_DOUBLE)
const TAB_Bxy_R = zeros(Float64, NV_DOUBLE, NV_DOUBLE)
const TAB_Byy_R = zeros(Float64, NV_DOUBLE, NV_DOUBLE)

##################################################
# Compute FFT of B
##################################################
function compute_B_FFT!()

    for i = (-NV+1):(NV-1)
        ux = i * DV
        for j = (-NV+1):(NV-1)
            uy = j * DV
            
            Bxx, Bxy, Byx, Byy = get_kernel_B(ux, uy)
            
            ind_i = mod(i, 2*NV) + 1
            ind_j = mod(j, 2*NV) + 1
            
            TAB_Bxx_R[ind_i, ind_j] = Bxx
            TAB_Bxy_R[ind_i, ind_j] = Bxy
            TAB_Byy_R[ind_i, ind_j] = Byy
        end
    end

    mul!(TAB_Bxx_K, PLAN_FFT_R_TO_K, TAB_Bxx_R)
    mul!(TAB_Bxy_K, PLAN_FFT_R_TO_K, TAB_Bxy_R)
    mul!(TAB_Byy_K, PLAN_FFT_R_TO_K, TAB_Byy_R)
end

# Initialising the Fourier collision kernel component, once and for all
compute_B_FFT!()
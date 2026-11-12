##################################################
# Compute Final Source term of S3 with Loop Fusion
# Computes (1/N) * \p_v2 . (\p_v3 . FLUX_S3)
##################################################

##################################################
# Thread-Local arrays for Zero-Allocation
##################################################
const TAB_DIV3_X_PARALLEL = [Matrix{Float64}(undef, NV2, NV2) for tid = 1:Threads.nthreads()+1]
const TAB_DIV3_Y_PARALLEL = [Matrix{Float64}(undef, NV2, NV2) for tid = 1:Threads.nthreads()+1]
#
const TAB_FLUX_XX_PARALLEL = [Vector{Float64}(undef, NV2) for tid = 1:Threads.nthreads()+1]
const TAB_FLUX_XY_PARALLEL = [Vector{Float64}(undef, NV2) for tid = 1:Threads.nthreads()+1]
const TAB_FLUX_YY_PARALLEL = [Vector{Float64}(undef, NV2) for tid = 1:Threads.nthreads()+1]
##################################################

function COMPUTE_TAB_SOURCE_S3!()
    fill!(TAB_SOURCE_S3, 0.0)
    
    precompute_S2_B_integral!()  # Compute int dv S2(v1, v) B(v, v2)
    
    # The external loop is l1 because derivatives \p_v2 and \p_v3 don't act on it.
    Threads.@threads :static for l1 in 1:NV2
        tid = Threads.threadid()
        # Local 2D buffers for the inner divergence field: \p_v3 . FLUX
        # Matrix is indexed as [l2, l3].
        local_div3_x  = TAB_DIV3_X_PARALLEL[tid]
        local_div3_y  = TAB_DIV3_Y_PARALLEL[tid]
        
        # Local 1D arrays for the Flux along v3 (for a fixed v1)        
        local_flux_xx = TAB_FLUX_XX_PARALLEL[tid]
        local_flux_xy = TAB_FLUX_XY_PARALLEL[tid]
        local_flux_yy = TAB_FLUX_YY_PARALLEL[tid]
        
        # Compute inner divergence wrt v3 for all (l2, l3)
        for l2 in 1:NV2
            # Read the variables depending only on l1 and l2 
            F_l2 = TAB_F[l2]
            int_F_Bxx = TAB_INT_Bxx_F[l2]
            int_F_Bxy = TAB_INT_Bxy_F[l2]
            int_F_Byy = TAB_INT_Byy_F[l2]
            
            S2_l2_l1 = TAB_S2[l2, l1]
            
            int_S2_Bxx = TAB_INT_S2_Bxx[l1, l2]
            int_S2_Bxy = TAB_INT_S2_Bxy[l1, l2]
            int_S2_Byy = TAB_INT_S2_Byy[l1, l2]
            
            tab_Bxx = TAB_TAB_Bxx[l2]
            tab_Bxy = TAB_TAB_Bxy[l2]
            tab_Byy = TAB_TAB_Byy[l2]
            
            # Fill the 1D Flux vectors along l3
            for l3 in 1:NV2
                #####
                F_l3 = TAB_F[l3]
                #
                Bxx_l2_l3 = tab_Bxx[l3]
                Bxy_l2_l3 = tab_Bxy[l3] 
                Byy_l2_l3 = tab_Byy[l3]
                #
                S2_l3_l1 = TAB_S2[l3, l1] 
                #
                dirac_l2_l3 = delta_cell(l2, l3)
                #
                local_flux_xx[l3] = - (S2_l2_l1 * F_l3 + S2_l3_l1 * F_l2) * Bxx_l2_l3 + 
                                     dirac_l2_l3 * ( S2_l2_l1 * int_F_Bxx + F_l2 * int_S2_Bxx )
                
                local_flux_xy[l3] = - (S2_l2_l1 * F_l3 + S2_l3_l1 * F_l2) * Bxy_l2_l3 + 
                                     dirac_l2_l3 * ( S2_l2_l1 * int_F_Bxy + F_l2 * int_S2_Bxy )
                
                local_flux_yy[l3] = - (S2_l2_l1 * F_l3 + S2_l3_l1 * F_l2) * Byy_l2_l3 + 
                                     dirac_l2_l3 * ( S2_l2_l1 * int_F_Byy + F_l2 * int_S2_Byy )
            end
            
            # Compute divergence wrt v3
            for l3 in 1:NV2
                ####
                dFxx = value_div_1D_field(local_flux_xx, 1, l3)
                dFxy = value_div_1D_field(local_flux_xy, 2, l3)
                local_div3_x[l2, l3] = dFxx + dFxy
                
                dFyx = value_div_1D_field(local_flux_xy, 1, l3)
                dFyy = value_div_1D_field(local_flux_yy, 2, l3)
                local_div3_y[l2, l3] = dFyx + dFyy
            end
        end
        
        # Compute outer divergence wrt v2
        for l3 in 1:NV2
            # Create a 1D view along l2 for a fixed l3
            div3_x_1D = @view local_div3_x[:, l3]
            div3_y_1D = @view local_div3_y[:, l3]
            
            for l2 in 1:NV2
                # Divergence wrt v2
                ddFxx = value_div_1D_field(div3_x_1D, 1, l2)
                ddFyy = value_div_1D_field(div3_y_1D, 2, l2)
                
                # Store the result in the final array
                # ATTENTION here ther is no factor 1/2
                TAB_SOURCE_S3[l1, l2, l3] = (ddFxx + ddFyy)
            end
        end
        
    end
end

##################################################
# Final Source term of S3 with the circular permutations
# ATTENTION, here we don't use the fact that S_3 is completely symmetric
##################################################
function compute_Source_S3_TOTAL!()

    # We compute the term with the order (l1, l2, l3)
    COMPUTE_TAB_SOURCE_S3!()
    
    # Compute the circular permutations
    Threads.@threads for l3 in 1:NV2
        for l2 in 1:NV2
            for l1 in 1:NV2
                # The first permutation (l1, l2, l3)
                term_123 = TAB_SOURCE_S3[l1, l2, l3]
                # Permutation SOURCE_S3(v2, v3, v1) 
                term_231 = TAB_SOURCE_S3[l2, l3, l1]
                # Permutation SOURCE_S3(v3, v1, v2) 
                term_312 = TAB_SOURCE_S3[l3, l1, l2]
    
                # Add all contributions
                TAB_dS3_dt[l1, l2, l3] += term_123 + term_231 + term_312
            end
        end
    end
end

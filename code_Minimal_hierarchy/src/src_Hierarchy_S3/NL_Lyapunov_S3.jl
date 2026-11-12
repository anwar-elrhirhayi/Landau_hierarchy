##################################################
# Compute Nonlinear S3 terms (S4 -> S2 S2)
##################################################

##################################################
# Thread-Local arrays for Zero-Allocation (NL Term)
##################################################
const TAB_NL_FLUX_X_PARALLEL = [Vector{Float64}(undef, NV2) for tid = 1:Threads.nthreads()+1]
const TAB_NL_FLUX_Y_PARALLEL = [Vector{Float64}(undef, NV2) for tid = 1:Threads.nthreads()+1]
##################################################

##################################################
# Compute the Flux inside \p_v1 and then the divergence
# of it
##################################################
function COMPUTE_TAB_NL_S3!()
    fill!(TAB_NL_S3, 0.0)
    
    # Computation of 1D and 2D integrals
    TAB_INTEGRALS_NL!()

    # Main loop over v2 and v3
    Threads.@threads :static for l2 in 1:NV2
        ####
        tid = Threads.threadid()
        # Allocate local buffers for the current thread.
        local_flux_x = TAB_NL_FLUX_X_PARALLEL[tid]
        local_flux_y = TAB_NL_FLUX_Y_PARALLEL[tid]
        
        for l3 in 1:NV2
            
            # We fill the flux with respect to v1 for the pair (l2, l3)
            for l1 in 1:NV2
                ##### Part (v2, v3)
                I2_xx_2 = TAB_INT_I2_xx[l1, l2]
                I2_xy_2 = TAB_INT_I2_xy[l1, l2]
                I2_yy_2 = TAB_INT_I2_yy[l1, l2]
                
                I3_x_2 = TAB_INT_I3_x[l1, l2]
                I3_y_2 = TAB_INT_I3_y[l1, l2]
                
                dS2_13_x = TAB_dS2_dx[l1, l3]
                dS2_13_y = TAB_dS2_dy[l1, l3]
                S2_13 = TAB_S2[l1, l3]
                
                Flux_23_x = (I2_xx_2 * dS2_13_x + I2_xy_2 * dS2_13_y) - S2_13 * I3_x_2
                Flux_23_y = (I2_xy_2 * dS2_13_x + I2_yy_2 * dS2_13_y) - S2_13 * I3_y_2
            
                ##### Part (v3, v2)
                I2_xx_3 = TAB_INT_I2_xx[l1, l3]
                I2_xy_3 = TAB_INT_I2_xy[l1, l3]
                I2_yy_3 = TAB_INT_I2_yy[l1, l3]
                
                I3_x_3 = TAB_INT_I3_x[l1, l3]
                I3_y_3 = TAB_INT_I3_y[l1, l3]
                
                dS2_12_x = TAB_dS2_dx[l1, l2]
                dS2_12_y = TAB_dS2_dy[l1, l2]
                S2_12 = TAB_S2[l1, l2]
                
                Flux_32_x = (I2_xx_3 * dS2_12_x + I2_xy_3 * dS2_12_y) - S2_12 * I3_x_3
                Flux_32_y = (I2_xy_3 * dS2_12_x + I2_yy_3 * dS2_12_y) - S2_12 * I3_y_3
                
                # Store in the local 1D vector
                local_flux_x[l1] = Flux_23_x + Flux_32_x
                local_flux_y[l1] = Flux_23_y + Flux_32_y
            end
            
            # Compute the divergence of this 1D vector
            for l1 in 1:NV2
                # Divergence 
                div_Flux_x = value_div_1D_field(local_flux_x, 1, l1)
                div_Flux_y = value_div_1D_field(local_flux_y, 2, l1)
                # 
                Term_23 = div_Flux_x + div_Flux_y
                # Add all contribution
                # ATTENTION to the factor 1/2
                TAB_NL_S3[l1, l2, l3] =  0.5 * Term_23 
            end
        end
    end
end

##################################################
# Apply circular permutations for total NL contribution
##################################################
function compute_NL_S3_TOTAL!()    
    ####
    COMPUTE_TAB_NL_S3!()
    ####
    Threads.@threads for l3 in 1:NV2
        for l2 in 1:NV2
            for l1 in 1:NV2
                # First term order (1, 2, 3)
                term_123 = TAB_NL_S3[l1, l2, l3]
                # Permutation (2, 3, 1)
                term_231 = TAB_NL_S3[l2, l3, l1]
                # Permutation (3, 1, 2)
                term_312 = TAB_NL_S3[l3, l1, l2]
                # Add all contributions
                TAB_dS3_dt[l1, l2, l3] += term_123 + term_231 + term_312
            end
        end
    end
end





##################################################
# OLD CODE
##################################################

##################################################
# Compute the Flux inside \p_v1
##################################################
#=
function TAB_FLUX_NL_S3!()
    fill!(TAB_FLUX_NL_x, 0.0)
    fill!(TAB_FLUX_NL_y, 0.0)
    
    Threads.@threads for l1 in 1:NV2
        for l2 in 1:NV2
            for l3 in 1:NV2
                ##### Part (v2, v3)
                # Tensor I2(v1, v2)
                I2_xx_2 = TAB_INT_I2_xx[l1, l2]
                I2_xy_2 = TAB_INT_I2_xy[l1, l2]
                I2_yy_2 = TAB_INT_I2_yy[l1, l2]
                # Vector I3(v1, v2)
                I3_x_2 = TAB_INT_I3_x[l1, l2]
                I3_y_2 = TAB_INT_I3_y[l1, l2]
                # \p_v1 S2(v1, v3)
                dS2_13_x = TAB_dS2_dx[l1, l3]
                dS2_13_y = TAB_dS2_dy[l1, l3]
                #
                S2_13 = TAB_S2[l1, l3]
                #
                Flux_23_x = (I2_xx_2 * dS2_13_x + I2_xy_2 * dS2_13_y) - S2_13 * I3_x_2
                Flux_23_y = (I2_xy_2 * dS2_13_x + I2_yy_2 * dS2_13_y) - S2_13 * I3_y_2
            
                ##### Part (v3, v2), using the symmetry (2, 3)
                I2_xx_3 = TAB_INT_I2_xx[l1, l3]
                I2_xy_3 = TAB_INT_I2_xy[l1, l3]
                I2_yy_3 = TAB_INT_I2_yy[l1, l3]
                #    
                I3_x_3 = TAB_INT_I3_x[l1, l3]
                I3_y_3 = TAB_INT_I3_y[l1, l3]
                #
                dS2_12_x = TAB_dS2_dx[l1, l2]
                dS2_12_y = TAB_dS2_dy[l1, l2]
                #
                S2_12 = TAB_S2[l1, l2]
                #
                Flux_32_x = (I2_xx_3 * dS2_12_x + I2_xy_3 * dS2_12_y) - S2_12 * I3_x_3
                Flux_32_y = (I2_xy_3 * dS2_12_x + I2_yy_3 * dS2_12_y) - S2_12 * I3_y_3
                # Sum the two contributions
                TAB_FLUX_NL_x[l1, l2, l3] = Flux_23_x + Flux_32_x
                TAB_FLUX_NL_y[l1, l2, l3] = Flux_23_y + Flux_32_y
            end
        end
    end
end
=#
##################################################
# Compute Divergences 
##################################################
#=
function COMPUTE_TAB_NL_S3!()
    fill!(TAB_NL_S3, 0.0)
    fill!(TAB_DIV_v1_FLUX_NL, 0.0)
    fill!(TAB_DIV_v1_I1, 0.0)
    
    TAB_INTEGRALS_NL!()
    TAB_FLUX_NL_S3!()
    
    # Divergence of the 1D field I1(v1)
    Threads.@threads for l1 in 1:NV2
        div_I1_x = value_div_1D_field(TAB_INT_I1_x, 1, l1)
        div_I1_y = value_div_1D_field(TAB_INT_I1_y, 2, l1)
        TAB_DIV_v1_I1[l1] = div_I1_x + div_I1_y
    end
    
    # Divergence of the 3D Flux Field \p_v1 . FLUX
    Threads.@threads for l1 in 1:NV2
        for l2 in 1:NV2
            for l3 in 1:NV2
                div_Flux_x = value_div_3D_field_v1(TAB_FLUX_NL_x, 1, l1, l2, l3)
                div_Flux_y = value_div_3D_field_v1(TAB_FLUX_NL_y, 2, l1, l2, l3)
                TAB_DIV_v1_FLUX_NL[l1, l2, l3] = div_Flux_x + div_Flux_y
            end
        end
    end
    
    Threads.@threads for l1 in 1:NV2
        for l2 in 1:NV2
            for l3 in 1:NV2
                S2_23 = TAB_S2[l2, l3]
                Term_1 = S2_23 * TAB_DIV_v1_I1[l1]
                Term_23 = TAB_DIV_v1_FLUX_NL[l1, l2, l3]
                #
                TAB_NL_S3[l1, l2, l3] = 0.5 * (Term_1 + Term_23)
            end
        end
    end
end
=#


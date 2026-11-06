##################################################
##################################################
# Lyapunov S_2^(1) and S_3^(2) 
##################################################
##################################################

##################################################
# Compute the RHS of Landau eq. p_t F(t) = RHS[F(t)]
##################################################
function compute_dF_dt_Hierarchy_S3!()
    ##### Computing the gradients
    TAB_dF!() 
    ##### Computing the Landau flux p_t F = div[Flux]
    TAB_FLUX_LANDAU!()    
    ##### Computing the divergence of the flux, i.e. dF/dt
    TAB_dF_dt!() 
end

##################################################
# Compute the RHS of Lyapunov eq. 
# p_t S_2^(1)(t) = RHS[F(t), S_2^(1)(t)]
# for given F(t) and S_2^(1)(t)
##################################################
function symmetrize_Lv_contribution_S2!()
    # Since S2(v1,v2) = S2(v2,v1), L_v2 S2 is the transpose of L_v1 S2.
    Threads.@threads :dynamic for l2 in 1:NV2
        @inbounds @simd for l1 in 1:(l2-1)
            Lv_sym = TAB_dS2_dt[l1, l2] + TAB_dS2_dt[l2, l1]
            TAB_dS2_dt[l1, l2] = Lv_sym
            TAB_dS2_dt[l2, l1] = Lv_sym
        end
        @inbounds TAB_dS2_dt[l2, l2] *= 2
    end
end

##################################################
# ATTENTION -- here we use the symmetry property of S2
##################################################
function compute_dS2_dt_Hierarchy_S3!()
    #####
    fill!(TAB_dS2_dt, 0.0)
    # Compute L_v1 contribution
    compute_Lv_contribution_S2!(TAB_F, TAB_dF_dx, TAB_dF_dy, 0)
    # Obtain L_v2 from L_v1 by symmetry
    symmetrize_Lv_contribution_S2!()
    # Source term
    compute_Source_S2!()
end

##################################################
# Compute the RHS of Lyapunov eq. for S3
# p_t S_3^(1)(t) = RHS[F(t), S_3^(2)(t)]
# for given F(t) and S_3^(2)(t)
##################################################
function complete_Lv_S3_from_symmetry!()
    # Save L_v1 before overwriting TAB_dS3_dt in TAB_NL_S3
    # ATTENTION -- here we use TAB_NL_S3 as a provisory container.
    # This operation is safe because  TAB_NL_S3 is filled with 0.0 juat becore computing
    # the non linear term.
    copyto!(TAB_NL_S3, TAB_dS3_dt)
    Threads.@threads :static for l3 in 1:NV2
        #####
        @inbounds for l2 in 1:NV2
            @simd for l1 in 1:NV2
                #####
                TAB_dS3_dt[l1, l2, l3] =
                    TAB_NL_S3[l1, l2, l3] +    # L_v1
                    TAB_NL_S3[l2, l1, l3] +    # L_v2
                    TAB_NL_S3[l3, l1, l2]      # L_v3
                #####
            end
        end
    end
end

function compute_dS3_dt_Hierarchy_S3!()
    #####
    fill!(TAB_dS3_dt, 0.0)
    # Compute grad S2 (needed for the NL term)
    TAB_dS2!()
    # Compute L_v1 once, then obtain L_v2 and L_v3 by symmetry
    compute_Lv_contribution_S3!(TAB_F, TAB_dF_dx, TAB_dF_dy, 0)
    # Symmetry for Lv_2 and Lv_3
    complete_Lv_S3_from_symmetry!()
    # Compute non linear contribution of S_4
    compute_NL_S3_TOTAL!()
    # Source term
    compute_Source_S3_TOTAL!()
end

##################################################
# Explicit Euler update
# Compute the drift F(t+dt) = F(t) + RHS[F(t)]*dt
##################################################
function drift_F_Hierarchy_S3!()
    for l in 1:NV2
        #####
        TAB_F[l] += DT * TAB_dF_dt[l] # Euler drift
        #####
    end
end

##################################################
# Explicit Euler update
# Compute the drift S2(t+dt) = S2(t) + RHS[F(t), S2(t)]*dt
##################################################
function drift_S2_Hierarchy_S3!()
    for l2 in 1:NV2,
        l1 in 1:NV2
        ####
        TAB_S2[l1, l2] += DT * TAB_dS2_dt[l1, l2] + DT * TAB_SOURCE_S2[l1, l2] # Euler drift
        ####
    end
end

##################################################
# Explicit Euler update
# Compute the drift S3(t+dt) = S3(t) + RHS[F(t), S2(t), S3(t)]*dt
##################################################
function drift_S3_Hierarchy_S3!()
    for l3 in 1:NV2,
        l2 in 1:NV2,
        l1 in 1:NV2
        ####
        TAB_S3[l1, l2, l3] += DT * TAB_dS3_dt[l1, l2, l3]
        ####
    end
end

##################################################
# Compute the RHS of Landau and Lyapunov eq. 
# for given F(t), S2(t) and S3(t)
# ATTENTION : RHS of S2 and S3 should be computed at 
# the same moment as F
##################################################
function compute_d_dt_Hierarchy_S3!()
    compute_dF_dt_Hierarchy_S3!()
    # Precompute integrals with F and B once and for all
    precompute_F_B_integrals!()
    #
    compute_dS2_dt_Hierarchy_S3!()
    compute_dS3_dt_Hierarchy_S3!()
end

##################################################
# Explicit Euler update
# Compute the drift 
# F(t+dt) = F(t) + RHS[F(t)]*dt
# S2(t+dt) = S2(t) + RHS[F(t), S2(t)]*dt
# S3(t+dt) = S3(t) + RHS[F(t), S2(t), S3(t)]*dt
##################################################
function drift_Hierarchy_S3!()
    drift_F_Hierarchy_S3!()
    drift_S2_Hierarchy_S3!()
    drift_S3_Hierarchy_S3!()
end

##################################################
# Integration of F, S2, S3 for one timestep.
##################################################
function integrate_dt_Hierarchy_S3!()
    compute_d_dt_Hierarchy_S3!()
    drift_Hierarchy_S3!()
end

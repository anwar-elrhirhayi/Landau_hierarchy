##################################################
##################################################
# Lyapunov S_2^(1) 
##################################################
##################################################

##################################################
# Compute the RHS of Landau eq. p_t F(t) = RHS[F(t)]
##################################################
function compute_dF_dt_Hierarchy_S2!()
    ##### Computing the gradients
    TAB_dF!()  # gradients of F
    TAB_dS2!() # gradients of S_2
    ##### Computing the flux p_t F = div[Flux]
    TAB_Flux_wto_correction!() 
    ##### The two functions here are used only to compute the Landau + 1/N S2 term
    # For the minimal hierarchy, this correction is not needed
    # TAB_Flux_correction_Hierarchy_S2!()  # 1/N-correction to the Landau-flux
    # TAB_Flux!()  # Computing the total flux
    ##### Computing the divergence of the flux, i.e. dF/dt
    TAB_dF_dt!() 
end

##################################################
# Explicit Euler update
# Compute the drift F(t+dt) = F(t) + RHS[F(t)]*dt
##################################################
function drift_F_Hierarchy_S2!()
    for l in 1:NV2
        #####
        TAB_F[l] += DT * TAB_dF_dt[l] # Euler drift
        #####
    end
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
function compute_dS2_dt_Hierarchy_S2!()
    #####
    fill!(TAB_dS2_dt, 0.0)
    # Precompute integrals with F in L_v1 once and for all
    precompute_F_integrals!()
    # Compute L_v1 contribution
    compute_Lv_contribution!(TAB_F, TAB_dF_dx, TAB_dF_dy, 0)
    # Obtain L_v2 from L_v1 by symmetry
    symmetrize_Lv_contribution_S2!()
    # Source term
    TAB_Q!()
end

function compute_dS2_dt_Hierarchy_S2_wto_symmetry!()
    #####
    fill!(TAB_dS2_dt, 0.0)
    # Precompute integrals with F in L_v1 once and for all
    precompute_F_integrals!()
    # Compute L_v1 contribution
    compute_Lv_contribution!(TAB_F, TAB_dF_dx, TAB_dF_dy, 0)
    # Compute L_v2 contribution
    compute_Lv_contribution!(TAB_F, TAB_dF_dx, TAB_dF_dy, 1)
    # Source term
    TAB_Q!()
end
##################################################
# Explicit Euler update
# Compute the drift S(t+dt) = S(t) + RHS[F(t), S(t)]*dt
##################################################
function drift_S2_Hierarchy_S2!()
    for l2  in 1:NV2,
        l1 in 1:NV2
        ####
        # ATTENTION : The factor 1/2 is important here
        # ATTENTION : No factor 1/N in front of TAB_Q
        TAB_S2[l1, l2] += 1/2 * DT * TAB_dS2_dt[l1, l2] + DT * TAB_Q[l1, l2] # Euler drift
        ####
    end
    ####
end

##################################################
# Compute the RHS of Landau and Lyapunov eq. 
# for given F(t) and S2(t)
# ATTENTION : RHS of S_2^(2) should be computed at 
# the same moment as F
##################################################
function compute_d_dt_Hierarchy_S2!()
    compute_dF_dt_Hierarchy_S2!()    
    compute_dS2_dt_Hierarchy_S2!()
end

##################################################
# Explicit Euler update
# Compute the drift 
# F(t+dt) = F(t) + RHS[F(t)]*dt
# S(t+dt) = S(t) + RHS[F(t), S(t)]*dt
##################################################
function drift_Hierarchy_S2!()
    drift_F_Hierarchy_S2!()
    drift_S2_Hierarchy_S2!()
end

##################################################
# Integration of S2 for one timestep.
##################################################
function integrate_dt_Hierarchy_S2!()
    compute_d_dt_Hierarchy_S2!()
    drift_Hierarchy_S2!()
end

##################################################
# Variance for S2_SF
##################################################
function compute_var_Tx_SF()
    #####
    sum_v2_F = 0.0
    sum_v4_F = 0.0
    #####
    for l in 1:NV2
        i = TAB_i[l]
        vx = TAB_Vx[i]
        f_val = TAB_F[l] 
        
        sum_v2_F += (vx^2) * f_val
        sum_v4_F += (vx^4) * f_val
    end
    #####
    integral_v2 = sum_v2_F * DV^2
    integral_v4 = sum_v4_F * DV^2
    
    res = integral_v4 - integral_v2^2
    
    # ATTENTION, the MASS_PART = 1/NPART factor comes here from the
    # fact that here we simulated S_2^{(1)} and VAR_Tx = < Tx, S_2 Tx'> = 1/NPART < Tx, S_2^(1) Tx'>
    res *= MASS_PART
    #####
    return res
end

function compute_var_Ty_SF()
    #####
    sum_vy2_F = 0.0
    sum_vy4_F = 0.0
    #####
    for l in 1:NV2
        j = TAB_j[l]
        vy = TAB_Vy[j]
        
        f_val = TAB_F[l] 
        
        sum_vy2_F += (vy^2) * f_val
        sum_vy4_F += (vy^4) * f_val
    end
    #####
    integral_vy2 = sum_vy2_F * DV^2
    integral_vy4 = sum_vy4_F * DV^2
    
    res = integral_vy4 - integral_vy2^2
    
    res *= MASS_PART
    #####
    return res
end

##################################################
# Gradients of F
# Finite differences: central in interior, one-sided at boundaries
# Results stored in TAB_dF_dx, TAB_dF_dy
##################################################
function TAB_dF!()
    compute_1D_field_gradients!(TAB_F,TAB_dF_dx,TAB_dF_dy) 
end

##################################################
# Divergence of Flux
# p_vx Flux_x + p_vy Flux_y (finite differences)
# Result stored in TAB_dF_dt
##################################################
function TAB_dF_dt!()
    #####
    # Fill TAB_dF_dt with div Flux
    compute_1D_field_div!(TAB_Fluxx, TAB_Fluxy, TAB_dF_dt)
    #####
end
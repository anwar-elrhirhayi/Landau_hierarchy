##################################################
# Gradients of S2
# Results stored in TAB_dS2_dx, TAB_dS2_dy
##################################################
function TAB_dS2!()
    compute_2D_field_gradients_v1(TAB_S2, TAB_dS2_dx, TAB_dS2_dy)  #  Along v1
end

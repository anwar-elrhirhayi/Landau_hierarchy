##################################################
# Gradients of S2
# Results stored in TAB_dS2_dx1, TAB_dS2_dy1, TAB_dS2_dx2, TAB_dS2_dy2
##################################################
function TAB_dS2!()
    compute_2D_field_gradients_v1(TAB_S2, TAB_dS2_dx1, TAB_dS2_dy1)  #  Along v1
    compute_2D_field_gradients_v2(TAB_S2, TAB_dS2_dx2, TAB_dS2_dy2)  #  Along v2
end

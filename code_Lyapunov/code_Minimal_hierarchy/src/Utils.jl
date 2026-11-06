##################################################
##################################################
# Utils for all Lyapunov hierarchy
##################################################
##################################################

##################################################
# Given (i,j) in the 2D tab --> l in the 1D tab
##################################################
function get_index_1D(i, j)
    return (i-1) * NV + j
end

##################################################
# Finite difference of a 1D field (1D array) on index v = (i, j) = l 
# If direction = 1 --> finite difference along vx direction
# If direction = 2 --> finite difference along vy direction
# Finite differences: central in interior, one-sided at boundaries
# Return a scalar TAB_DIFF_ij
# Result: value of the finite difference evaluated in v = l
##################################################

function compute_finite_diff_at_v(TAB_1D, l, direction)
    TAB_DIFF_ij = 0.0
    i, j = TAB_i[l], TAB_j[l]
    
    if direction == 1 # derivative on vx = i (direction == 1)
        ####
        if i == 1 
            TAB_DIFF_ij = (-3 * TAB_1D[l] + 4 * TAB_1D[l + NV] - TAB_1D[l + 2*NV]) / (2 * DV)
        elseif i == NV
            TAB_DIFF_ij = (3 * TAB_1D[l] - 4 * TAB_1D[l - NV] + TAB_1D[l - 2*NV]) / (2 * DV)
        else
            TAB_DIFF_ij = (TAB_1D[l + NV] - TAB_1D[l - NV]) / (2 * DV)
        end
        
    elseif direction == 2 # derivative on vy = j (direction == 2)
        ####
        if j == 1 
            TAB_DIFF_ij = (-3 * TAB_1D[l] + 4 * TAB_1D[l + 1] - TAB_1D[l + 2]) / (2 * DV)
        elseif j == NV
            TAB_DIFF_ij = (3 * TAB_1D[l] - 4 * TAB_1D[l - 1] + TAB_1D[l - 2]) / (2 * DV)
        else
            TAB_DIFF_ij = (TAB_1D[l + 1] - TAB_1D[l - 1]) / (2 * DV) 
        end
    end
    
    return TAB_DIFF_ij
end

#=
function compute_finite_diff_at_v(TAB_1D, l, direction)
    TAB_DIFF_ij = 0.0
    i, j = TAB_i[l], TAB_j[l]
    if     direction == 1 # derivative on vx = i (direction == 1)
        ####
        if i == 1 
            TAB_DIFF_ij = (TAB_1D[l + NV] - TAB_1D[l]) / (2 * DV)      # Boundary V_x = -VMAX
        elseif i == NV
            TAB_DIFF_ij = (TAB_1D[l]      - TAB_1D[l - NV]) / (2 * DV) # Boundary V_x =  VMAX
        else
            TAB_DIFF_ij = (TAB_1D[l + NV] - TAB_1D[l - NV]) / (2 * DV)
        end
        
    elseif direction == 2 # derivative on vy = j (direction == 2)
        ####
        if j == 1 
            TAB_DIFF_ij = (TAB_1D[l + 1] - TAB_1D[l]) /  (2 * DV)      # Boundary V_y = -VMAX
        elseif j == NV
            TAB_DIFF_ij = (TAB_1D[l]     - TAB_1D[l - 1]) /  (2 * DV)  # Boundary V_y =  VMAX
        else
            TAB_DIFF_ij = (TAB_1D[l + 1] - TAB_1D[l - 1]) / (2 * DV) 
        end
    end
    return TAB_DIFF_ij
end
=#
##################################################
# Divergence of a 1D field flux on v = (i,j) = l
# Result: value of the divergence evaluated in v = l
##################################################

function value_div_1D_field(Flux, direction, l)
    divF = 0.0
    if direction == 1
        ##### # derivative on vx = i (direction == 1)
        i = TAB_i[l]
        if i == 1 
            divF = (3 * Flux[l] + Flux[l+NV]) / (2 * DV)
        elseif i == 2
            divF = (-4 * Flux[l-NV] + Flux[l+NV]) / (2 * DV)
        elseif i == 3
            divF = (Flux[l-2*NV] - Flux[l-NV] + Flux[l+NV]) / (2 * DV)
        elseif i == NV
            divF = (-Flux[l-NV] - 3 * Flux[l]) / (2 * DV)
        elseif i == NV - 1
            divF = (-Flux[l-NV] + 4 * Flux[l+NV]) / (2 * DV)
        elseif i == NV - 2
            divF = (-Flux[l-NV] + Flux[l+NV] - Flux[l+2*NV]) / (2 * DV)
        else
            divF = (Flux[l+NV] - Flux[l-NV]) / (2 * DV)
        end
        #####
    elseif direction == 2
        ##### # derivative on vy = j (direction == 2)
        j = TAB_j[l]
        if j == 1 
            divF = (3 * Flux[l] + Flux[l+1]) / (2 * DV)
        elseif j == 2
            divF = (-4 * Flux[l-1] + Flux[l+1]) / (2 * DV)
        elseif j == 3
            divF = (Flux[l-2] - Flux[l-1] + Flux[l+1]) / (2 * DV)
        elseif j == NV
            divF = (-Flux[l-1] - 3 * Flux[l]) / (2 * DV)
        elseif j == NV - 1
            divF = (-Flux[l-1] + 4 * Flux[l+1]) / (2 * DV)
        elseif j == NV - 2
            divF = (-Flux[l-1] + Flux[l+1] - Flux[l+2]) / (2 * DV)
        else
            divF = (Flux[l+1] - Flux[l-1]) / (2 * DV)
        end
        #####
    end
    return divF
end
#=
function value_div_1D_field(Flux, direction, l)
    divF = 0.0
    if direction == 1
        #####
        i = TAB_i[l]
        if i == 1
            divF = (Flux[l+NV] + Flux[l]) / (2 * DV)
        elseif i == NV
            divF = - (Flux[l] + Flux[l-NV]) / (2 * DV)
        else
            divF = (Flux[l+NV] - Flux[l-NV]) / (2 * DV)
        end
        #####
    elseif direction == 2
        #####
        j = TAB_j[l]
        if j == 1
            divF = (Flux[l+1] + Flux[l]) / (2 * DV)
        elseif j == NV
            divF = - (Flux[l] + Flux[l-1]) / (2 * DV)
        else
            divF = (Flux[l+1] - Flux[l-1]) / (2 * DV)
        end
        #####
    end
    return divF
end
=#

##################################################
# Compute gradient field (dGdx, dGdy) of a 1D field G
# by central differences 
# results in dGdx[l], dGdy[l]
##################################################
function compute_1D_field_gradients!(G, dGdx, dGdy)
    for l in 1:NV2
        # derivative on vx = i (direction == 1)
        direction = 1
        dGdx[l] = compute_finite_diff_at_v(G, l, direction)
        # derivative on vy = j (direction == 2)
        direction = 2
        dGdy[l] = compute_finite_diff_at_v(G, l, direction)
    end
end

##################################################
# Gradient of a 2D field with respect to v1, evaluate on
# a given point v1, v2 = l1, l2
##################################################
function value_gradient_2D_field_v1(field, direction, l1, l2)
    ####
    tab_1D = @view field[:, l2]
    ####
    dfield = compute_finite_diff_at_v(tab_1D, l1, direction)
    return dfield
end

##################################################
# Gradient of a 2D field with respect to v2, evaluate on
# a given point v1, v2 = l1, l2
####################################
function value_gradient_2D_field_v2(field, direction, l1, l2)
    ####
    tab_1D = @view field[l1, :]
    ####
    dfield = compute_finite_diff_at_v(tab_1D, l2, direction)
    return dfield
end

##################################################
# Given some input 2D field (2 variables array) tab_in,
# and output 2 1D fields (1D arrays), tab_grad_x, tab_grad_y,
# it computes the gradient fields on x an y.
##################################################

# Gradients along the first variable v1
function compute_2D_field_gradients_v1(tab_in, tab_grad_x1, tab_grad_y1)
    for l2 = 1:NV2,
        l1 = 1:NV2
        #####
        tab_grad_x1[l1, l2] = value_gradient_2D_field_v1(tab_in, 1, l1, l2)
        tab_grad_y1[l1, l2] = value_gradient_2D_field_v1(tab_in, 2, l1, l2)
        #####
    end
    return tab_grad_x1, tab_grad_y1
end
#
# Gradients along the second variable v2
function compute_2D_field_gradients_v2(tab_in, tab_grad_x2, tab_grad_y2)
    for l2 = 1:NV2,
        l1 = 1:NV2
        #####
        tab_grad_x2[l1, l2] = value_gradient_2D_field_v2(tab_in, 1, l1, l2)
        tab_grad_y2[l1, l2] = value_gradient_2D_field_v2(tab_in, 2, l1, l2)
        #####
    end
    return tab_grad_x2, tab_grad_y2
end

##################################################
##################################################

##################################################
# Divergence of a 2D field with respect to v1, evaluated
# in a given point (v1, v2) = (l1, l2)
##################################################
function value_div_2D_field_v1(Flux_2D, direction, l1, l2)
    ####
    Flux_1D = @view Flux_2D[:, l2]
    ####
    divFlux = value_div_1D_field(Flux_1D, direction, l1)
    return divFlux
end

##################################################
# Divergence of a 2D field with respect to v2,
# evaluated in a given point (v1, v2) = (l1, l2)
##################################################
function value_div_2D_field_v2(Flux_2D, direction, l1, l2)
    ####
    Flux_1D = @view Flux_2D[l1, :]
    ####
    divFlux = value_div_1D_field(Flux_1D, direction, l2)
    return divFlux
end

##################################################
# Compute 1D field of div of a 1D field
##################################################
function compute_1D_field_div!(tab_in_x, tab_in_y, tab_div)
    for l = 1:NV2
        #####
        div_x = value_div_1D_field(tab_in_x, 1, l)
        div_y = value_div_1D_field(tab_in_y, 2, l)
        #####
        tab_div[l] = div_x + div_y
    end
end


##################################################
# Extracts a 1D field S_field(v) from the 2D field TAB_S2(v, v')
# depending on the targeted variable
# indx = 0 --> v' fixed
# indx = 1 --> v  fixed
# Return S_field evaluated on (k, l):
# If indx = 0, return S_field(.) = S(., v'), v' = k fixed
# If indx = 1, return S_field(.) = S(v, .),  v  = k fixed
##################################################
function slice_2D_field(indx, k, TAB_S2_SLICE)
    if     indx == 0 # we fix v' = k
        for l in 1:NV2
            TAB_S2_SLICE[l] = TAB_S2[l, k]
        end
    elseif indx == 1 # we fix v  = k
        for l in 1:NV2
            TAB_S2_SLICE[l] = TAB_S2[k, l]
        end
    end
end  

##################################################
# Slicing a 3D field TAB_S3 into a 1D field TAB_S3_SLICE
# depending on the targeted variable (indx = 0, 1 or 2)
##################################################
function slice_3D_field(indx, l2, l3, TAB_S3_SLICE)
    if indx == 0
        # We fix the 2nd and 3rd variables to l2 and l3
        # S3_slice(v) = S3(v, v2, v3)
        @inbounds for l in 1:NV2
            TAB_S3_SLICE[l] = TAB_S3[l, l2, l3]
        end
    elseif indx == 1
        # We fix the 1st and 3rd variables to l2 and l3
        # S3_slice(v) = S3(v3, v, v2)
        @inbounds for l in 1:NV2
            TAB_S3_SLICE[l] = TAB_S3[l3, l, l2]
        end
    elseif indx == 2
        # We fix the 1st and 2nd variables to l2 and l3
        # S3_slice(v) = S3(v2, v3, v)
        @inbounds for l in 1:NV2
            TAB_S3_SLICE[l] = TAB_S3[l2, l3, l]
        end
    end
end


##################################################
# Dirac(v1 - v2) where
# v1 = l1 and v2 = l2
##################################################
function delta_cell(l1, l2)
    if l1 == l2
        return 1.0 / (DV2) 
    else
        return 0.0
    end
end

##################################################
# Regularized Dirac delta 
# to allow stable finite-difference derivatives
##################################################
# function delta_cell_smooth(l1, l2)
#     i1, j1 = TAB_i[l1], TAB_j[l1]
#     i2, j2 = TAB_i[l2], TAB_j[l2]
#     v1_x = TAB_Vx[i1]
#     v1_y = TAB_Vy[j1]
    
#     v2_x = TAB_Vx[i2]
#     v2_y = TAB_Vy[j2]
#     sigma = 1.5 * sqrt(DV2) 
#     dist_sq = (v1_x - v2_x)^2 + (v1_y - v2_y)^2
#     norm_factor = 1.0 / (pi * sigma^2)
    
#     return norm_factor * exp(-dist_sq / (sigma^2))
# end

##################################################
# Divergence of a 3D field with respect to v1, v2 or v3
# evaluated in a given point (v1, v2, v3) = (l1, l2, l3)
##################################################
function value_div_3D_field_v3(Flux_3D, direction, l1, l2, l3)
    ####
    Flux_1D = @view Flux_3D[l1, l2, :]
    ####
    divFlux = value_div_1D_field(Flux_1D, direction, l3)
    return divFlux
end
#
function value_div_3D_field_v2(Flux_3D, direction, l1, l2, l3)
    ####
    Flux_1D = @view Flux_3D[l1, :, l3]
    ####
    divFlux = value_div_1D_field(Flux_1D, direction, l2)
    return divFlux
end
#
function value_div_3D_field_v1(Flux_3D, direction, l1, l2, l3)
    ####
    Flux_1D = @view Flux_3D[:, l2, l3]
    ####
    divFlux = value_div_1D_field(Flux_1D, direction, l1)
    return divFlux
end
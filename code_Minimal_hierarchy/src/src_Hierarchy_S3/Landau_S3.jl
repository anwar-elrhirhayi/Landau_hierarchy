
##################################################
# Landau equation in 2D + correction term
# Discretization of velocity space with double integral
# Explicit Euler scheme for time integration
##################################################
# Compute the flux, Flux(v1):
# p_t F = div Flux(v1)
# Flux(v1) = p_v1 int dv2 B(v1, v2) [- p_v2 + p_v1] [F(v1) F(v2) + 1/N S_2^(1)(v1,v2)]
# Implementation: double loop over v1 and v2
# Results stored in TAB_Fluxx, TAB_Fluxy
##################################################

##################################################
# Flux without correction
##################################################
function TAB_FLUX_LANDAU!()
    ##### Initialising the arrays
    fill!(TAB_Fluxx, 0.0)
    fill!(TAB_Fluxy, 0.0)
    
    #####
    Threads.@threads for l1 in 1:NV2
        # Read the TAB
        @inbounds begin
            Fv1    = TAB_F[l1]
            dF_dx1 = TAB_dF_dx[l1]
            dF_dy1 = TAB_dF_dy[l1]
        end
        
        # Initialize the accumulated collision terms in x and y directions for particle i
        acc_x = 0.0 
        acc_y = 0.0 
        
        #####
        # Loop over v2 = l2
        #####
        @inbounds @simd for l2 in 1:NV2
            #####
            # ATTENTION, here we used the symmetry of B[l1, l2]
            # We read [l2, l1] instead of [l1, l2] 
            # to read TAB_B sequentially (column by column)
            Bxx = TAB_Bxx[l2, l1]
            Bxy = TAB_Bxy[l2, l1]
            Byx = Bxy
            Byy = TAB_Byy[l2, l1]
            
            #####
            dF_dx2 = TAB_dF_dx[l2]
            dF_dy2 = TAB_dF_dy[l2]
            Fv2     = TAB_F[l2]
            
            #####
            # Compute integrand
            # pref = B(v1, v2) [- F(v1) p_v2 F + F(v2) p_v1 F]
            #####
            pref_x = @muladd(-dF_dx2 * Fv1 + Fv2 * dF_dx1)
            pref_y = @muladd(-dF_dy2 * Fv1 + Fv2 * dF_dy1)
            
            #####
            term_x = @muladd(Bxx * pref_x + Bxy * pref_y)
            term_y = @muladd(Byx * pref_x + Byy * pref_y)
            
            #####
            acc_x += term_x
            acc_y += term_y
        end
        
        ##### Multiplying by the volume element
        # ATTENTION to the factor 1/2 here
        acc_x *= DV^2 * 1/2
        acc_y *= DV^2 * 1/2
        
        ##### Filling in the array
        TAB_Fluxx[l1] = acc_x
        TAB_Fluxy[l1] = acc_y
    end
end

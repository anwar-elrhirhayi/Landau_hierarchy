##################################################
# Landau equation in 2D
# Discretization of velocity space with double integral
# Explicit Euler scheme for time integration
##################################################
# Compute the flux, Flux(v):
# p_t F = div Flux(v)
# Flux(v) = p_v int dv' B(v, v') [- F(v) p_v'F + F(v') p_vF]
# Implementation: double loop over v and v'
# Results stored in TAB_Fluxx, TAB_Fluxy
##################################################
function TAB_Flux!()
    ##### Initialising the arrays
    fill!(TAB_Fluxx, 0.0)
    fill!(TAB_Fluxy, 0.0)
    
    #####
    Threads.@threads for l in 1:NV2
        # Read the TAB
        @inbounds begin
            Fv    = TAB_F[l]
            dF_dx = TAB_dF_dx[l]
            dF_dy = TAB_dF_dy[l]
        end
        
        # Initialize the accumulated collision terms in x and y directions for particle i
        acc_x = 0.0 
        acc_y = 0.0 
        
        #####
        # Loop over v' = lp
        #####
        @inbounds @simd for lp in 1:NV2
            #####
            # ATTENTION, here we used the symmetry of B[l, lp]
            # We read [lp, l] instead of [l, lp] 
            # to read TAB_B sequentially (column by column)
            Bxx = TAB_Bxx[lp, l]
            Bxy = TAB_Bxy[lp, l]
            Byx = Bxy
            Byy = TAB_Byy[lp, l]
            
            #####
            dF_dxp = TAB_dF_dx[lp]
            dF_dyp = TAB_dF_dy[lp]
            Fp     = TAB_F[lp]
            
            #####
            # Compute integrand
            # pref = B(v, v') [- F(v) p_v'F + F(v') p_vF]
            #####
            pref_x = @muladd(-dF_dxp * Fv + Fp * dF_dx)
            pref_y = @muladd(-dF_dyp * Fv + Fp * dF_dy)
            
            #####
            term_x = @muladd(Bxx * pref_x + Bxy * pref_y)
            term_y = @muladd(Byx * pref_x + Byy * pref_y)
            
            #####
            acc_x += term_x
            acc_y += term_y
        end
        
        ##### Multiplying by the volume element
        acc_x *= DV^2 * 1/2
        acc_y *= DV^2 * 1/2
        
        ##### Filling in the array
        @inbounds begin
            TAB_Fluxx[l] = acc_x
            TAB_Fluxy[l] = acc_y
        end
    end
end

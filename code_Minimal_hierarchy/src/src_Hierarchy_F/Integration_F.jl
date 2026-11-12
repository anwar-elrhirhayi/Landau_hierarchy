##################################################
##################################################
# Landau eq
##################################################
##################################################

##################################################
# Compute the RHS of Landau eq. p_t F(t) = RHS[F(t)]
##################################################
function compute_dF_dt!()
    ##### Computing the gradients
    TAB_dF!()       
    ##### Computing the flux p_t F = div[Flux]
    TAB_Flux!()               
    ##### Computing the divergence of the flux, i.e. dF/dt
    TAB_dF_dt!() 
end

##################################################
# Explicit Euler update
# Compute the drift F(t+dt) = F(t) + RHS[F(t)]*dt
##################################################
function drift_F!()
    for l in 1:NV2
        #####
        TAB_F[l] += DT * TAB_dF_dt[l] # Euler drift
        #####
    end
end

##################################################
# Integration of F for one timestep.
##################################################
function integrate_dt_hierarchy_F!()
    compute_dF_dt!() # Computing the rates of changes
    drift_F!()       # Pushing forward in time
end

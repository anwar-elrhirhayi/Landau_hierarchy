##################################################
# Initial conditions for S2^(1)
##################################################
function init_S2!()
    for l1 in 1:NV2
        F_l1 = TAB_F[l1]
        #####
        for l2 in 1:NV2
            # ATTENTION, no factor 1/N 
            TAB_S2[l1, l2] = delta_cell(l1, l2) * F_l1 - TAB_F[l2] * F_l1
        end
    end
end

##################################################
# Initialisising on F --> S
##################################################
function init_Hierarchy_S2!(Tx, Ty)
    init_F!(Tx, Ty)
    init_S2!()
end

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
# Initial conditions for S3^(2)
##################################################
function init_S3!()
    for l1 in 1:NV2
        F_l1 = TAB_F[l1]
        #####
        for l2 in 1:NV2
            F_l2 = TAB_F[l2]
            d_12 = delta_cell(l1, l2)
            for l3 = 1:NV2
                #####
                F_l3 = TAB_F[l3]
                d_23 = delta_cell(l2, l3)
                d_13 = delta_cell(l1, l3)
                # ATTENTION, no factor 1/N^2 
                TAB_S3[l1, l2, l3] = (
                    d_12 * d_23 * F_l1 
                  - F_l1 * F_l3 * d_12
                  - F_l1 * F_l2 * d_23
                  - F_l2 * F_l3 * d_13 
                  + 2.0 * F_l1 * F_l2 * F_l3
                )
            end
        end
    end
end

##################################################
# Initialisising on F --> S
##################################################
function init_Hierarchy_S3!(Tx, Ty)
    init_F!(Tx, Ty)
    #####
    init_S2!()
    init_S3!()
end
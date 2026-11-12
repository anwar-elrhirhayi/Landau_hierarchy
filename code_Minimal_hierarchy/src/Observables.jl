##################################################
# Total mass
##################################################
function compute_mass()
    #####
    res = 0.0
    #####
    for l in 1:NV2
        #####
        res += TAB_F[l]
        #####
    end
    #####
    res *= DV^(2)
    #####
    return res
end
##################################################
# Total momentum
##################################################
function compute_momentum()
    #####
    res_x = 0.0
    res_y = 0.0
    #####
    for l in 1:NV2
        #####
        i, j = TAB_i[l], TAB_j[l]
        #####
        vx = TAB_Vx[i]
        vy = TAB_Vy[j]
        vF = TAB_F[l]
        res_x += vx * vF
        res_y += vy * vF
        #####
    end
    #####
    res_x *= DV^(2)
    res_y *= DV^(2)
    #####
    return res_x, res_y
end
#
function compute_Ptot()
    Ptot1, Ptot2 = compute_momentum()
    #####
    return sqrt(Ptot1^2 + Ptot2^2)
end
##################################################
# Total energy
##################################################
function compute_energy()
    #####
    res = 0.0
    #####
    for l in 1:NV2
        #####
        i, j = TAB_i[l], TAB_j[l]
        #####
        vx = TAB_Vx[i]
        vy = TAB_Vy[j]
        vF = TAB_F[l]
        vE = vx^(2) + vy^(2)
        res += vE * vF
        #####
    end
    #####
    res *= DV^(2)
    #####
    return res
end

##################################################
# Var(MASS) = \int dv1 dv2 S(v1, v2)
##################################################
function compute_var_mass()
    #####
    res = 0.0
    #####
    for l1 in 1:NV2
        for l2 in 1:NV2
            #####
            res += TAB_S2[l1, l2]
            #####
        end
    end
    #####
    res *= DV^(4) * MASS_PART
    #####
    return res
end

function compute_err_var_mass()
    return abs((compute_var_mass() - TAB_VAR_M[1])/(TAB_VAR_M[1]))
end
##################################################
# Var(Ekin) = \int dv1 dv2 v1^2 S(v1, v2) v2^2 
##################################################
function compute_var_energy()
    #####
    res = 0.0
    #####
    for l1 in 1:NV2
        #####
        i1, j1 = TAB_i[l1], TAB_j[l1]
        #####
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]
        for l2 in 1:NV2
            i2, j2 = TAB_i[l2], TAB_j[l2]
            #####
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            res += (vx1^2 + vy1^2) * TAB_S2[l1, l2] * (vx2^2 + vy2^2)
            #####
        end
    end
    #####
    res *= DV^(4) * MASS_PART
    #####
    return res
end

function compute_err_var_Etot()
    return abs((compute_var_energy() - TAB_VAR_E[1])/(TAB_VAR_E[1]))
end
##################################################
# Var(P_x) = \int dv1 dv2 vx1 S(v1, v2) vx2
#          = 1/N \int dv1 dv2 vx1 S^(1)(v1, v2) vx2
##################################################
function compute_var_Px()
    #####
    res = 0.0
    #####
    for l1 in 1:NV2
        i1 = TAB_i[l1]
        vx1 = TAB_Vx[i1]
        for l2 in 1:NV2
            i2 = TAB_i[l2]
            #####
            vx2 = TAB_Vx[i2]
            res += vx1 * vx2 * TAB_S2[l1, l2]
            #####
        end
    end
    #####
    # ATTENTION to the MASS_PART = 1/NPART factor 
    res *= DV^(4) * MASS_PART
    #####
    return res
end

##################################################
# Var(P_y) = \int dv1 dv2 vy1 S(v1, v2) vy2 
##################################################
function compute_var_Py()
    #####
    res = 0.0
    #####
    for l1 in 1:NV2
        j1 = TAB_j[l1]
        vy1 = TAB_Vy[j1]
        for l2 in 1:NV2
            j2 = TAB_j[l2]
            #####
            vy2 = TAB_Vy[j2]
            res += vy1 * vy2 * TAB_S2[l1, l2]
            #####
        end
    end
    #####
    res *= DV^(4) * MASS_PART
    #####
    return res
end

##################################################
# Var(Ptot) = \int dv1 dv2 Ptot(v1) S(v1, v2) Ptot(v2) 
##################################################
function compute_var_Ptot()
    #####
    res = 0.0
    #####
    for l1 in 1:NV2
        #####
        i1, j1 = TAB_i[l1], TAB_j[l1]
        #####
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]
        for l2 in 1:NV2
            i2, j2 = TAB_i[l2], TAB_j[l2]
            #####
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            res += sqrt(vx1^2 + vy1^2) * TAB_S2[l1, l2] * sqrt(vx2^2 + vy2^2)
            #####
        end
    end
    #####
    res *= DV^(4) * MASS_PART
    #####
    return res
end

###############################################
# Error in the total energy and momentum
# ATTENTION, this is a RELATIVE error
###############################################
function compute_err_Etot()
    return abs((compute_energy() - TAB_E[1])/(TAB_E[1]))
end
#
function compute_err_Ptot()
    Ptot1, Ptot2 = compute_momentum()
    #####
    norm_DeltaPtot = sqrt((Ptot1 - TAB_Px[1])^(2) +
                          (Ptot2 - TAB_Py[1])^(2))
    #####
    norm_Ptotinit = sqrt(TAB_Px[1]^(2) + TAB_Py[1]^(2))
    #####
    err_Ptot = norm_DeltaPtot / norm_Ptotinit
    #####
    return err_Ptot
end

##################################################
# M3(MASS) = \int dv1 dv2 dv3 S(v1, v2, v3)
##################################################
function compute_M3_mass()
    #####
    res = 0.0
    #####
    @inbounds for l1 in 1:NV2
        for l2 in 1:NV2            
            for l3 in 1:NV2
                #####
                res += TAB_S3[l1, l2, l3]
                #####
            end
        end
    end
    #####
    # Triple integral in 2D, Vol = DV^6
    # M3 scales as 1/N^2 = MASS_PART^2
    res *= DV^(6) * MASS_PART^2
    #####
    return res
end

function compute_err_M3_mass()
    return abs((compute_M3_mass() - TAB_M3_M[1])/(TAB_M3_M[1]))
end
##################################################
# M3(Ptot) = \int dv1 dv2 dv3 Ptot(v1) Ptot(v2) Ptot(v3) S_3(v1, v2, v3)
##################################################
function compute_M3_Px_Py()
    #####
    res_x = 0.0
    res_y = 0.0
    #####
    for l1 in 1:NV2
        #####
        i1, j1 = TAB_i[l1], TAB_j[l1]
        #
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]
        #
        for l2 in 1:NV2
            #####
            i2, j2 = TAB_i[l2], TAB_j[l2]
            #
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            #
            for l3 in 1:NV2
                #####
                i3, j3 = TAB_i[l3], TAB_j[l3]
                #
                vx3 = TAB_Vx[i3]
                vy3 = TAB_Vy[j3]
                #
                res_x += vx1 * vx2 * vx3 * TAB_S3[l1, l2, l3] 
                res_y += vy1 * vy2 * vy3 * TAB_S3[l1, l2, l3] 
            end
        end
    end
    #####
    res_x *= DV^(6) * MASS_PART^2
    res_y *= DV^(6) * MASS_PART^2
    #####
    return res_x, res_y
end

function compute_M3_Ptot()
    #####
    res = 0.0
    #####
    for l1 in 1:NV2
        #####
        i1, j1 = TAB_i[l1], TAB_j[l1]
        #
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]
        Ptot1 = sqrt(vx1^2 + vy1^2)
        #
        for l2 in 1:NV2
            #####
            i2, j2 = TAB_i[l2], TAB_j[l2]
            #
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            Ptot2 = sqrt(vx2^2 + vy2^2)
            for l3 in 1:NV2
                #####
                i3, j3 = TAB_i[l3], TAB_j[l3]
                #
                vx3 = TAB_Vx[i3]
                vy3 = TAB_Vy[j3]
                Ptot3 = sqrt(vx3^2 + vy3^2)
                #
                res += Ptot1 * Ptot2 * Ptot3 * TAB_S3[l1, l2, l3] 
            end
        end
    end
    #####
    res *= DV^(6) * MASS_PART^2
    #####
    return res
end

##################################################
# M3(Ekin)
##################################################

function compute_M3_energy()
    #####
    res = 0.0
    #####
    for l1 in 1:NV2
        #####
        i1, j1 = TAB_i[l1], TAB_j[l1]
        #
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]
        E1 = vx1^2 + vy1^2
        #
        for l2 in 1:NV2
            #####
            i2, j2 = TAB_i[l2], TAB_j[l2]
            #
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            E2 = vx2^2 + vy2^2
            for l3 in 1:NV2
                #####
                i3, j3 = TAB_i[l3], TAB_j[l3]
                #
                vx3 = TAB_Vx[i3]
                vy3 = TAB_Vy[j3]
                E3 = vx3^2 + vy3^2
                #
                res += E1 * E2 * E3 * TAB_S3[l1, l2, l3] 
            end
        end
    end
    #####
    res *= DV^(6) * MASS_PART^2
    #####
    return res
end

function compute_err_M3_Etot()
    return abs((compute_M3_energy() - TAB_M3_E[1])/(TAB_M3_E[1]))
end
##################################################
# Observable: temperatures
# Tx = 1/MASS_TOT * \sum m vx_i^2  
##################################################
##################################################
function compute_Tx()
    #####
    res = 0.0
    #####
    for l in 1:NV2
        #####
        i = TAB_i[l]
        #####
        vx = TAB_Vx[i]
        res += vx^2 * TAB_F[l]
        #####
    end
    #####
    res *= DV^(2)
    #####
    return res
end
##################################################
# Ty = int vy^2 F dv
##################################################
function compute_Ty()
    #####
    res = 0.0
    #####
    for l in 1:NV2
        #####
        j = TAB_j[l]
        #####
        vy = TAB_Vy[j]
        res += vy^2 * TAB_F[l]
        #####
    end
    #####
    res *= DV^(2)
    #####
    return res
end

##################################################
# Var(T_x) = \int dv1 dv2 vx1^2 S(v1, v2) vx2^2
#          = 1/N \int dv1 dv2 vx1^2 S^(1)(v1, v2) vx2^2
##################################################
function compute_var_Tx()
    #####
    res = 0.0
    #####
    for l1 in 1:NV2
        i1 = TAB_i[l1]
        vx1 = TAB_Vx[i1]
        for l2 in 1:NV2
            i2 = TAB_i[l2]
            #####
            vx2 = TAB_Vx[i2]
            res += vx1^2 * vx2^2 * TAB_S2[l1, l2]
            #####
        end
    end
    #####
    # ATTENTION, the MASS_PART = 1/NPART factor comes here froms the
    # fact that here we simulated S_2^{(1)} and VAR_Tx = < Tx, S_2 Tx'> = 1/NPART < Tx, S_2^(1) Tx'>
    res *= DV^(4) * MASS_PART
    #####
    return res
end

##################################################
# Var(T_y) = \int dv1 dv2 vy1^2 S(v1, v2) vy2^2 
##################################################
function compute_var_Ty()
    #####
    res = 0.0
    #####
    for l1 in 1:NV2
        j1 = TAB_j[l1]
        vy1 = TAB_Vy[j1]
        for l2 in 1:NV2
            j2 = TAB_j[l2]
            #####
            vy2 = TAB_Vy[j2]
            res += vy1^2 * vy2^2 * TAB_S2[l1, l2]
            #####
        end
    end
    #####
    res *= DV^(4) * MASS_PART
    #####
    return res
end

##################################################
# M3(T_x) = \int dv1 dv2 dv3 vx1^2 vx2^2 vx3^2 S3(v1, v2, v3)
##################################################
function compute_M3_Tx()
    #####
    res = 0.0
    #####
    @inbounds for l1 in 1:NV2
        i1 = TAB_i[l1]
        vx1_sq = TAB_Vx[i1]^2
        #
        for l2 in 1:NV2
            i2 = TAB_i[l2]
            #####
            vx2_sq = TAB_Vx[i2]^2
            vx12_sq = vx1_sq * vx2_sq 
            
            for l3 in 1:NV2
                i3 = TAB_i[l3]
                #####
                vx3_sq = TAB_Vx[i3]^2
                res += vx12_sq * vx3_sq * TAB_S3[l1, l2, l3]
                #####
            end
        end
    end
    #####
    # Triple integral in 2D, Vol = DV^6
    # M3 scales as 1/N^2 = MASS_PART^2
    res *= DV^(6) * MASS_PART^2
    #####
    return res
end

##################################################
# M3(T_y) = \int dv1 dv2 dv3 vy1^2 vy2^2 vy3^2 S3(v1, v2, v3)
##################################################
function compute_M3_Ty()
    #####
    res = 0.0
    #####
    @inbounds for l1 in 1:NV2
        j1 = TAB_j[l1]
        vy1_sq = TAB_Vy[j1]^2
        #
        for l2 in 1:NV2
            j2 = TAB_j[l2]
            #####
            vy2_sq = TAB_Vy[j2]^2
            vy12_sq = vy1_sq * vy2_sq 
            #####
            for l3 in 1:NV2
                j3 = TAB_j[l3]
                #####
                vy3_sq = TAB_Vy[j3]^2
                res += vy12_sq * vy3_sq * TAB_S3[l1, l2, l3]
                #####
            end
        end
    end
    #####
    res *= DV^(6) * MASS_PART^2
    #####
    # var_Ty = compute_var_Ty()
    # standardized_skew_Ty = res / (var_Ty^1.5)
    return res
end


##################################################
# Compute the centred 4th momentum using the closure on S_4^(2)
# MOM4(T_x) = \int S4(1, 2, 3, 4) T_x(1) ... T_x(4)
# at order 1/N^2
##################################################

function compute_delta_initial_moment4_Tx()
    #####
    mu2 = 0.0
    mu4 = 0.0
    mu6 = 0.0
    mu8 = 0.0
    #####
    for l in 1:NV2
        #####
        i = TAB_i[l]
        #
        vx = TAB_Vx[i]
        F0 = TAB_F[l]
        #
        vx2 = vx^2
        mu2 += vx2 * F0
        mu4 += vx2^2 * F0
        mu6 += vx2^3 * F0
        mu8 += vx2^4 * F0
        #####
    end
    #####
    mu2 *= DV^(2)
    mu4 *= DV^(2)
    mu6 *= DV^(2)
    mu8 *= DV^(2)
    #####
    # Connected part of S_4
    # S_4^(3)(t=0) = 1/N^3 * kappa_4
    kappa4_part = mu8 - 4 * mu6 * mu2 - 3 * mu4^2 + 12 * mu4 * mu2^2 - 6 * mu2^4
    #####
    res = kappa4_part * MASS_PART^(3)
    #####
    return res
end

function compute_delta_initial_moment4_Ty()
    #####
    mu2 = 0.0
    mu4 = 0.0
    mu6 = 0.0
    mu8 = 0.0
    #####
    for l in 1:NV2
        #####
        j = TAB_j[l]
        #
        vy = TAB_Vy[j]
        F0 = TAB_F[l]
        #
        vy2 = vy^2
        mu2 += vy2 * F0
        mu4 += vy2^2 * F0
        mu6 += vy2^3 * F0
        mu8 += vy2^4 * F0
        #####
    end
    #####
    mu2 *= DV^(2)
    mu4 *= DV^(2)
    mu6 *= DV^(2)
    mu8 *= DV^(2)
    #####
    # Connected part of S_4
    # S_4^(3)(t=0) = 1/N^3 * kappa_4
    kappa4_part = mu8 - 4 * mu6 * mu2 - 3 * mu4^2 + 12 * mu4 * mu2^2 - 6 * mu2^4
    #####
    res = kappa4_part * MASS_PART^(3)
    #####
    return res
end

##################################################
# Cubic Observable the velocities
# i.e., S_xxx = 1 / MTOT * \sum_i^N m_i * vx_i * vx_i * vx_i
#       S_xxy = 1 / MTOT * \sum_i^N m_i * vx_i * vx_i * vy_i
#       S_xyy = 1 / MTOT * \sum_i^N m_i * vx_i * vy_i * vy_i
#       S_yyy = 1 / MTOT * \sum_i^N m_i * vy_i * vy_i * vy_i
###############################################
function compute_cubic()
    #####
    S_xxx, S_xxy, S_xyy, S_yyy = 0.0, 0.0, 0.0, 0.0
    #####
    for l in 1:NV2
        #####
        i, j = TAB_i[l], TAB_j[l]
        #####
        vx = TAB_Vx[i]
        vy = TAB_Vy[j]
        f_val = TAB_F[l]
        #####
        S_xxx += vx^3 * f_val
        S_xxy += vx^2 * vy * f_val
        S_xyy += vx * vy^2 * f_val
        S_yyy += vy^3 * f_val
        #####
    end
    #####
    S_xxx *= DV^(2)
    S_xxy *= DV^(2)
    S_xyy *= DV^(2)
    S_yyy *= DV^(2)
    #####
    return S_xxx, S_xxy, S_xyy, S_yyy
end

function compute_var_cubic()
    #####
    var_Sxxx, var_Sxxy, var_Sxyy, var_Syyy = 0.0, 0.0, 0.0, 0.0
    #####
    @inbounds for l1 in 1:NV2
        i1, j1 = TAB_i[l1], TAB_j[l1]      
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]

        for l2 in 1:NV2
            i2, j2 = TAB_i[l2], TAB_j[l2]  
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            #####
            s2_val = TAB_S2[l1, l2]
            
            var_Sxxx += vx1^3 * vx2^3 * s2_val
            var_Sxxy += vx1^2 * vy1 * s2_val * vx2^2 * vy2
            var_Sxyy += vx1 * vy1^2 * s2_val * vx2 * vy2^2
            var_Syyy += vy1^3 * s2_val * vy2^3
            #####
        end
    end
    #####
    # ATTENTION, the MASS_PART = 1/NPART factor comes here from the
    # fact that here we simulated S_2^{(1)} and VAR = 1/NPART \int S S_2 S' >
    factor = DV^(4) * MASS_PART
    
    var_Sxxx *= factor
    var_Sxxy *= factor
    var_Sxyy *= factor
    var_Syyy *= factor
    #####
    return var_Sxxx, var_Sxxy, var_Sxyy, var_Syyy
end

function compute_M3_cubic()
    #####
    m3_Sxxx, m3_Sxxy, m3_Sxyy, m3_Syyy = 0.0, 0.0, 0.0, 0.0
    #####
    @inbounds for l1 in 1:NV2
        i1, j1 = TAB_i[l1], TAB_j[l1]  
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]
        
        vx1_3 = vx1^3
        vx1_2_vy1 = vx1^2 * vy1
        vx1_vy1_2 = vx1 * vy1^2
        vy1_3 = vy1^3
        
        for l2 in 1:NV2
            i2, j2 = TAB_i[l2], TAB_j[l2]
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            
            vx12_3 = vx1_3 * vx2^3
            vx12_2_vy = vx1_2_vy1 * (vx2^2 * vy2)
            vx12_vy_2 = vx1_vy1_2 * (vx2 * vy2^2)
            vy12_3 = vy1_3 * vy2^3
            
            for l3 in 1:NV2
                i3, j3 = TAB_i[l3], TAB_j[l3]
                vx3 = TAB_Vx[i3]
                vy3 = TAB_Vy[j3]
                #####
                s3_val = TAB_S3[l1, l2, l3]
                
                m3_Sxxx += vx12_3 * vx3^3 * s3_val
                m3_Sxxy += vx12_2_vy * (vx3^2 * vy3) * s3_val
                m3_Sxyy += vx12_vy_2 * (vx3 * vy3^2) * s3_val
                m3_Syyy += vy12_3 * vy3^3 * s3_val
                #####
            end
        end
    end
    #####
    # Triple integral in 2D, Vol = DV^6
    # M3 scales as 1/N^2 = MASS_PART^2
    factor = DV^(6) * MASS_PART^2
    
    m3_Sxxx *= factor
    m3_Sxxy *= factor
    m3_Sxyy *= factor
    m3_Syyy *= factor
    #####
    return m3_Sxxx, m3_Sxxy, m3_Sxyy, m3_Syyy
end

##################################################
# Non polynomial Observable the velocities
##################################################
function observable_non_polyn(vx, vy)
    if vx >= 0
        return 1
    else
        return 0
    end
end

function compute_avg_observable_non_polyn()
    Obs_val = 0
    #####
    for l in 1:NV2
        #####
        i, j = TAB_i[l], TAB_j[l]
        #
        vx = TAB_Vx[i]
        vy = TAB_Vy[j]
        f_val = TAB_F[l]
        #
        Obs_val += f_val * observable_non_polyn(vx, vy) 
        ####
    end
    Obs_val *= DV^2 
    return Obs_val
end

function compute_var_observable_non_polyn()
    #####
    Obs_var = 0
    #####
    @inbounds for l1 in 1:NV2
        #####
        i1, j1 = TAB_i[l1], TAB_j[l1]      
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]
        #
        Obs_1 = observable_non_polyn(vx1, vy1)
        #####
        for l2 in 1:NV2
            i2, j2 = TAB_i[l2], TAB_j[l2]  
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            #
            Obs_2 = observable_non_polyn(vx2, vy2)
            #####
            s2_val = TAB_S2[l1, l2]
            Obs_var += Obs_1 * s2_val * Obs_2
        end
    end
    Obs_var *= DV^(4) * MASS_PART
    #####
    return Obs_var
end

function compute_M3_observable_non_polyn()
    #####
    m3_obs = 0.0
    #####
    @inbounds for l1 in 1:NV2
        #####
        i1, j1 = TAB_i[l1], TAB_j[l1]  
        #####
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]
        #
        obs_1 = observable_non_polyn(vx1, vy1)
        #####
        for l2 in 1:NV2
            #####
            i2, j2 = TAB_i[l2], TAB_j[l2]
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            #
            obs_2 = observable_non_polyn(vx2, vy2)
            for l3 in 1:NV2
                #####
                i3, j3 = TAB_i[l3], TAB_j[l3]
                vx3 = TAB_Vx[i3]
                vy3 = TAB_Vy[j3]
                #
                obs_3 = observable_non_polyn(vx3, vy3)
                #####
                m3_obs += TAB_S3[l1, l2, l3] * obs_1 * obs_2 * obs_3
                #####
            end
        end
    end
    #####
    # Triple integral in 2D, Vol = DV^6
    # M3 scales as 1/N^2 = MASS_PART^2
    factor = DV^(6) * MASS_PART^2
    #
    m3_obs *= factor
    #####
    return m3_obs
end

##################################################
# Quartic Observable the velocities
# i.e., Q_xxxx = 1 / MTOT * \sum_i^N m_i * vx_i * vx_i * vx_i * vx_i
#       Q_xxxy = 1 / MTOT * \sum_i^N m_i * vx_i * vx_i * vx_i * vy_i
#       Q_xxyy = 1 / MTOT * \sum_i^N m_i * vx_i * vx_i * vy_i * vy_i
#       Q_xyyy = 1 / MTOT * \sum_i^N m_i * vx_i * vy_i * vy_i * vy_i
#       Q_yyyy = 1 / MTOT * \sum_i^N m_i * vy_i * vy_i * vy_i * vy_i
##################################################
function compute_quartic()
    #####
    Q_xxxx, Q_xxxy, Q_xxyy, Q_xyyy, Q_yyyy = 0.0, 0.0, 0.0, 0.0, 0.0
    #####
    for l in 1:NV2
        #####
        i, j = TAB_i[l], TAB_j[l]
        #####
        vx = TAB_Vx[i]
        vy = TAB_Vy[j]
        f_val = TAB_F[l]
        #####
        Q_xxxx += vx^4 * f_val
        Q_xxxy += vx^3 * vy * f_val
        Q_xxyy += vx^2 * vy^2 * f_val
        Q_xyyy += vx * vy^3 * f_val
        Q_yyyy += vy^4 * f_val
        #####
    end
    #####
    Q_xxxx *= DV^(2)
    Q_xxxy *= DV^(2)
    Q_xxyy *= DV^(2)
    Q_xyyy *= DV^(2)
    Q_yyyy *= DV^(2)
    #####
    return Q_xxxx, Q_xxxy, Q_xxyy, Q_xyyy, Q_yyyy
end

function compute_var_quartic()
    #####
    var_Qxxxx, var_Qxxxy, var_Qxxyy, var_Qxyyy, var_Qyyyy = 0.0, 0.0, 0.0, 0.0, 0.0
    #####
    @inbounds for l1 in 1:NV2
        i1, j1 = TAB_i[l1], TAB_j[l1]
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]

        for l2 in 1:NV2
            i2, j2 = TAB_i[l2], TAB_j[l2]
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            #####
            s2_val = TAB_S2[l1, l2]

            var_Qxxxx += vx1^4 * vx2^4 * s2_val
            var_Qxxxy += vx1^3 * vy1 * s2_val * vx2^3 * vy2
            var_Qxxyy += vx1^2 * vy1^2 * s2_val * vx2^2 * vy2^2
            var_Qxyyy += vx1 * vy1^3 * s2_val * vx2 * vy2^3
            var_Qyyyy += vy1^4 * s2_val * vy2^4
            #####
        end
    end
    #####
    # ATTENTION, the MASS_PART = 1/NPART factor comes here from the
    # fact that here we simulated S_2^{(1)} and VAR = 1/NPART \int S S_2 S' >
    factor = DV^(4) * MASS_PART

    var_Qxxxx *= factor
    var_Qxxxy *= factor
    var_Qxxyy *= factor
    var_Qxyyy *= factor
    var_Qyyyy *= factor
    #####
    return var_Qxxxx, var_Qxxxy, var_Qxxyy, var_Qxyyy, var_Qyyyy
end

function compute_M3_quartic()
    #####
    m3_Qxxxx, m3_Qxxxy, m3_Qxxyy, m3_Qxyyy, m3_Qyyyy = 0.0, 0.0, 0.0, 0.0, 0.0
    #####
    @inbounds for l1 in 1:NV2
        i1, j1 = TAB_i[l1], TAB_j[l1]
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]

        vx1_4 = vx1^4
        vx1_3_vy1 = vx1^3 * vy1
        vx1_2_vy1_2 = vx1^2 * vy1^2
        vx1_vy1_3 = vx1 * vy1^3
        vy1_4 = vy1^4

        for l2 in 1:NV2
            i2, j2 = TAB_i[l2], TAB_j[l2]
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]

            vx12_4 = vx1_4 * vx2^4
            vx12_3_vy = vx1_3_vy1 * (vx2^3 * vy2)
            vx12_2_vy_2 = vx1_2_vy1_2 * (vx2^2 * vy2^2)
            vx12_vy_3 = vx1_vy1_3 * (vx2 * vy2^3)
            vy12_4 = vy1_4 * vy2^4

            for l3 in 1:NV2
                i3, j3 = TAB_i[l3], TAB_j[l3]
                vx3 = TAB_Vx[i3]
                vy3 = TAB_Vy[j3]
                #####
                s3_val = TAB_S3[l1, l2, l3]

                m3_Qxxxx += vx12_4 * vx3^4 * s3_val
                m3_Qxxxy += vx12_3_vy * (vx3^3 * vy3) * s3_val
                m3_Qxxyy += vx12_2_vy_2 * (vx3^2 * vy3^2) * s3_val
                m3_Qxyyy += vx12_vy_3 * (vx3 * vy3^3) * s3_val
                m3_Qyyyy += vy12_4 * vy3^4 * s3_val
                #####
            end
        end
    end
    #####
    # Triple integral in 2D, Vol = DV^6
    # M3 scales as 1/N^2 = MASS_PART^2
    factor = DV^(6) * MASS_PART^2

    m3_Qxxxx *= factor
    m3_Qxxxy *= factor
    m3_Qxxyy *= factor
    m3_Qxyyy *= factor
    m3_Qyyyy *= factor
    #####
    return m3_Qxxxx, m3_Qxxxy, m3_Qxxyy, m3_Qxyyy, m3_Qyyyy
end

##################################################
# Poly-sin-cos Observable the velocities
# i.e., Obs = 1 / MTOT * \sum_i^N m_i * vx_i * sin(vx_i) * vy_i^2 * sin(vy_i)
##################################################
function poly_sin_cos_obs(vx, vy)
    return vx * sin(vx) * vy^2 * sin(vy)
end

function compute_avg_poly_sin_cos()
    Obs_val = 0.0
    #####
    for l in 1:NV2
        #####
        i, j = TAB_i[l], TAB_j[l]
        #
        vx = TAB_Vx[i]
        vy = TAB_Vy[j]
        f_val = TAB_F[l]
        #
        Obs_val += f_val * poly_sin_cos_obs(vx, vy)
        ####
    end
    Obs_val *= DV^2
    return Obs_val
end

function compute_var_poly_sin_cos()
    #####
    Obs_var = 0.0
    #####
    @inbounds for l1 in 1:NV2
        #####
        i1, j1 = TAB_i[l1], TAB_j[l1]
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]
        #
        Obs_1 = poly_sin_cos_obs(vx1, vy1)
        #####
        for l2 in 1:NV2
            i2, j2 = TAB_i[l2], TAB_j[l2]
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            #
            Obs_2 = poly_sin_cos_obs(vx2, vy2)
            #####
            s2_val = TAB_S2[l1, l2]
            Obs_var += Obs_1 * s2_val * Obs_2
        end
    end
    Obs_var *= DV^(4) * MASS_PART
    #####
    return Obs_var
end

function compute_M3_poly_sin_cos()
    #####
    m3_obs = 0.0
    #####
    @inbounds for l1 in 1:NV2
        #####
        i1, j1 = TAB_i[l1], TAB_j[l1]
        #####
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]
        #
        obs_1 = poly_sin_cos_obs(vx1, vy1)
        #####
        for l2 in 1:NV2
            #####
            i2, j2 = TAB_i[l2], TAB_j[l2]
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            #
            obs_2 = poly_sin_cos_obs(vx2, vy2)
            for l3 in 1:NV2
                #####
                i3, j3 = TAB_i[l3], TAB_j[l3]
                vx3 = TAB_Vx[i3]
                vy3 = TAB_Vy[j3]
                #
                obs_3 = poly_sin_cos_obs(vx3, vy3)
                #####
                m3_obs += TAB_S3[l1, l2, l3] * obs_1 * obs_2 * obs_3
                #####
            end
        end
    end
    #####
    # Triple integral in 2D, Vol = DV^6
    # M3 scales as 1/N^2 = MASS_PART^2
    factor = DV^(6) * MASS_PART^2
    #
    m3_obs *= factor
    #####
    return m3_obs
end

##################################################
# Gaussian Observable the velocities
# i.e., Obs = 1 / MTOT * \sum_i^N m_i * exp(-(vx_i^2 + vy_i^2))
##################################################
function gaussian_obs(vx, vy)
    return exp(-(vx^2 + vy^2))
end

function compute_avg_gaussian()
    Obs_val = 0.0
    #####
    for l in 1:NV2
        #####
        i, j = TAB_i[l], TAB_j[l]
        #
        vx = TAB_Vx[i]
        vy = TAB_Vy[j]
        f_val = TAB_F[l]
        #
        Obs_val += f_val * gaussian_obs(vx, vy)
        ####
    end
    Obs_val *= DV^2
    return Obs_val
end

function compute_var_gaussian()
    #####
    Obs_var = 0.0
    #####
    @inbounds for l1 in 1:NV2
        #####
        i1, j1 = TAB_i[l1], TAB_j[l1]
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]
        #
        Obs_1 = gaussian_obs(vx1, vy1)
        #####
        for l2 in 1:NV2
            i2, j2 = TAB_i[l2], TAB_j[l2]
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            #
            Obs_2 = gaussian_obs(vx2, vy2)
            #####
            s2_val = TAB_S2[l1, l2]
            Obs_var += Obs_1 * s2_val * Obs_2
        end
    end
    Obs_var *= DV^(4) * MASS_PART
    #####
    return Obs_var
end

function compute_M3_gaussian()
    #####
    m3_obs = 0.0
    #####
    @inbounds for l1 in 1:NV2
        #####
        i1, j1 = TAB_i[l1], TAB_j[l1]
        #####
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]
        #
        obs_1 = gaussian_obs(vx1, vy1)
        #####
        for l2 in 1:NV2
            #####
            i2, j2 = TAB_i[l2], TAB_j[l2]
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            #
            obs_2 = gaussian_obs(vx2, vy2)
            for l3 in 1:NV2
                #####
                i3, j3 = TAB_i[l3], TAB_j[l3]
                vx3 = TAB_Vx[i3]
                vy3 = TAB_Vy[j3]
                #
                obs_3 = gaussian_obs(vx3, vy3)
                #####
                m3_obs += TAB_S3[l1, l2, l3] * obs_1 * obs_2 * obs_3
                #####
            end
        end
    end
    #####
    # Triple integral in 2D, Vol = DV^6
    # M3 scales as 1/N^2 = MASS_PART^2
    factor = DV^(6) * MASS_PART^2
    #
    m3_obs *= factor
    #####
    return m3_obs
end

##################################################
# Anisotropic channel pair: cos(vx) +/- cos(vy)
# O_plus  : isotropic test (symmetric under x vs y)
# O_minus : anisotropic test (antisymmetric, drives T_x != T_y relaxation)
##################################################
function cos_plus_obs(vx, vy)
    return cos(vx) + cos(vy)
end

function cos_minus_obs(vx, vy)
    return cos(vx) - cos(vy)
end

function compute_avg_cos_plus()
    Obs_val = 0.0
    #####
    for l in 1:NV2
        #####
        i, j = TAB_i[l], TAB_j[l]
        #
        vx = TAB_Vx[i]
        vy = TAB_Vy[j]
        f_val = TAB_F[l]
        #
        Obs_val += f_val * cos_plus_obs(vx, vy)
        ####
    end
    Obs_val *= DV^2
    return Obs_val
end

function compute_avg_cos_minus()
    Obs_val = 0.0
    #####
    for l in 1:NV2
        #####
        i, j = TAB_i[l], TAB_j[l]
        #
        vx = TAB_Vx[i]
        vy = TAB_Vy[j]
        f_val = TAB_F[l]
        #
        Obs_val += f_val * cos_minus_obs(vx, vy)
        ####
    end
    Obs_val *= DV^2
    return Obs_val
end

function compute_var_cos_plus()
    #####
    Obs_var = 0.0
    #####
    @inbounds for l1 in 1:NV2
        #####
        i1, j1 = TAB_i[l1], TAB_j[l1]
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]
        #
        Obs_1 = cos_plus_obs(vx1, vy1)
        #####
        for l2 in 1:NV2
            i2, j2 = TAB_i[l2], TAB_j[l2]
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            #
            Obs_2 = cos_plus_obs(vx2, vy2)
            #####
            s2_val = TAB_S2[l1, l2]
            Obs_var += Obs_1 * s2_val * Obs_2
        end
    end
    Obs_var *= DV^(4) * MASS_PART
    #####
    return Obs_var
end

function compute_var_cos_minus()
    #####
    Obs_var = 0.0
    #####
    @inbounds for l1 in 1:NV2
        #####
        i1, j1 = TAB_i[l1], TAB_j[l1]
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]
        #
        Obs_1 = cos_minus_obs(vx1, vy1)
        #####
        for l2 in 1:NV2
            i2, j2 = TAB_i[l2], TAB_j[l2]
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            #
            Obs_2 = cos_minus_obs(vx2, vy2)
            #####
            s2_val = TAB_S2[l1, l2]
            Obs_var += Obs_1 * s2_val * Obs_2
        end
    end
    Obs_var *= DV^(4) * MASS_PART
    #####
    return Obs_var
end

function compute_M3_cos_plus()
    #####
    m3_obs = 0.0
    #####
    @inbounds for l1 in 1:NV2
        #####
        i1, j1 = TAB_i[l1], TAB_j[l1]
        #####
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]
        #
        obs_1 = cos_plus_obs(vx1, vy1)
        #####
        for l2 in 1:NV2
            #####
            i2, j2 = TAB_i[l2], TAB_j[l2]
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            #
            obs_2 = cos_plus_obs(vx2, vy2)
            for l3 in 1:NV2
                #####
                i3, j3 = TAB_i[l3], TAB_j[l3]
                vx3 = TAB_Vx[i3]
                vy3 = TAB_Vy[j3]
                #
                obs_3 = cos_plus_obs(vx3, vy3)
                #####
                m3_obs += TAB_S3[l1, l2, l3] * obs_1 * obs_2 * obs_3
                #####
            end
        end
    end
    #####
    factor = DV^(6) * MASS_PART^2
    #
    m3_obs *= factor
    #####
    return m3_obs
end

function compute_M3_cos_minus()
    #####
    m3_obs = 0.0
    #####
    @inbounds for l1 in 1:NV2
        #####
        i1, j1 = TAB_i[l1], TAB_j[l1]
        #####
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]
        #
        obs_1 = cos_minus_obs(vx1, vy1)
        #####
        for l2 in 1:NV2
            #####
            i2, j2 = TAB_i[l2], TAB_j[l2]
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            #
            obs_2 = cos_minus_obs(vx2, vy2)
            for l3 in 1:NV2
                #####
                i3, j3 = TAB_i[l3], TAB_j[l3]
                vx3 = TAB_Vx[i3]
                vy3 = TAB_Vy[j3]
                #
                obs_3 = cos_minus_obs(vx3, vy3)
                #####
                m3_obs += TAB_S3[l1, l2, l3] * obs_1 * obs_2 * obs_3
                #####
            end
        end
    end
    #####
    factor = DV^(6) * MASS_PART^2
    #
    m3_obs *= factor
    #####
    return m3_obs
end

##################################################
# Observable O(v) = vx*sin(vx) * vy*sin(vy)
##################################################
function vsin_vsin_obs(vx, vy)
    return (vx * sin(vx)) * (vy * sin(vy))
end

function compute_avg_vsin_vsin()
    Obs_val = 0.0
    #####
    for l in 1:NV2
        #####
        i, j = TAB_i[l], TAB_j[l]
        #
        vx = TAB_Vx[i]
        vy = TAB_Vy[j]
        f_val = TAB_F[l]
        #
        Obs_val += f_val * vsin_vsin_obs(vx, vy)
        ####
    end
    Obs_val *= DV^2
    return Obs_val
end

function compute_var_vsin_vsin()
    #####
    Obs_var = 0.0
    #####
    @inbounds for l1 in 1:NV2
        #####
        i1, j1 = TAB_i[l1], TAB_j[l1]
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]
        #
        Obs_1 = vsin_vsin_obs(vx1, vy1)
        #####
        for l2 in 1:NV2
            i2, j2 = TAB_i[l2], TAB_j[l2]
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            #
            Obs_2 = vsin_vsin_obs(vx2, vy2)
            #####
            s2_val = TAB_S2[l1, l2]
            Obs_var += Obs_1 * s2_val * Obs_2
        end
    end
    Obs_var *= DV^(4) * MASS_PART
    #####
    return Obs_var
end

function compute_M3_vsin_vsin()
    #####
    m3_obs = 0.0
    #####
    @inbounds for l1 in 1:NV2
        #####
        i1, j1 = TAB_i[l1], TAB_j[l1]
        #####
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]
        #
        obs_1 = vsin_vsin_obs(vx1, vy1)
        #####
        for l2 in 1:NV2
            #####
            i2, j2 = TAB_i[l2], TAB_j[l2]
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            #
            obs_2 = vsin_vsin_obs(vx2, vy2)
            for l3 in 1:NV2
                #####
                i3, j3 = TAB_i[l3], TAB_j[l3]
                vx3 = TAB_Vx[i3]
                vy3 = TAB_Vy[j3]
                #
                obs_3 = vsin_vsin_obs(vx3, vy3)
                #####
                m3_obs += TAB_S3[l1, l2, l3] * obs_1 * obs_2 * obs_3
                #####
            end
        end
    end
    #####
    factor = DV^(6) * MASS_PART^2
    #
    m3_obs *= factor
    #####
    return m3_obs
end

##################################################
# Bounded odd observable in vx: O = tanh(vx)
##################################################
function tanh_vx_obs(vx)
    return tanh(vx)
end

function compute_avg_tanh_vx()
    Obs_val = 0.0
    #####
    for l in 1:NV2
        #####
        i = TAB_i[l]
        #
        vx = TAB_Vx[i]
        f_val = TAB_F[l]
        #
        Obs_val += f_val * tanh_vx_obs(vx)
        ####
    end
    Obs_val *= DV^2
    return Obs_val
end

function compute_var_tanh_vx()
    #####
    Obs_var = 0.0
    #####
    @inbounds for l1 in 1:NV2
        i1 = TAB_i[l1]
        vx1 = TAB_Vx[i1]
        #
        Obs_1 = tanh_vx_obs(vx1)
        #####
        for l2 in 1:NV2
            i2 = TAB_i[l2]
            #####
            vx2 = TAB_Vx[i2]
            Obs_2 = tanh_vx_obs(vx2)
            #####
            s2_val = TAB_S2[l1, l2]
            Obs_var += Obs_1 * s2_val * Obs_2
        end
    end
    Obs_var *= DV^(4) * MASS_PART
    #####
    return Obs_var
end

function compute_M3_tanh_vx()
    #####
    m3_obs = 0.0
    #####
    @inbounds for l1 in 1:NV2
        #####
        i1 = TAB_i[l1]
        #
        vx1 = TAB_Vx[i1]
        #
        obs_1 = tanh_vx_obs(vx1)
        #####
        for l2 in 1:NV2
            #####
            i2 = TAB_i[l2]
            vx2 = TAB_Vx[i2]
            #
            obs_2 = tanh_vx_obs(vx2)
            for l3 in 1:NV2
                #####
                i3 = TAB_i[l3]
                vx3 = TAB_Vx[i3]
                #
                obs_3 = tanh_vx_obs(vx3)
                #####
                m3_obs += TAB_S3[l1, l2, l3] * obs_1 * obs_2 * obs_3
                #####
            end
        end
    end
    #####
    factor = DV^(6) * MASS_PART^2
    #
    m3_obs *= factor
    #####
    return m3_obs
end

##################################################
# Bounded odd observable in vx 
# O = vx / sqrt(1 + vx^2)
##################################################
function sigmoid_vx_obs(vx)
    return vx / sqrt(1.0 + vx * vx)
end

function compute_avg_sigmoid_vx()
    Obs_val = 0.0
    #####
    for l in 1:NV2
        #####
        i = TAB_i[l]
        #
        vx = TAB_Vx[i]
        f_val = TAB_F[l]
        #
        Obs_val += f_val * sigmoid_vx_obs(vx)
        ####
    end
    Obs_val *= DV^2
    return Obs_val
end

function compute_var_sigmoid_vx()
    #####
    Obs_var = 0.0
    #####
    @inbounds for l1 in 1:NV2
        i1 = TAB_i[l1]
        vx1 = TAB_Vx[i1]
        #
        Obs_1 = sigmoid_vx_obs(vx1)
        #####
        for l2 in 1:NV2
            i2 = TAB_i[l2]
            #####
            vx2 = TAB_Vx[i2]
            Obs_2 = sigmoid_vx_obs(vx2)
            #####
            s2_val = TAB_S2[l1, l2]
            Obs_var += Obs_1 * s2_val * Obs_2
        end
    end
    Obs_var *= DV^(4) * MASS_PART
    #####
    return Obs_var
end

function compute_M3_sigmoid_vx()
    #####
    m3_obs = 0.0
    #####
    @inbounds for l1 in 1:NV2
        #####
        i1 = TAB_i[l1]
        #
        vx1 = TAB_Vx[i1]
        #
        obs_1 = sigmoid_vx_obs(vx1)
        #####
        for l2 in 1:NV2
            #####
            i2 = TAB_i[l2]
            vx2 = TAB_Vx[i2]
            #
            obs_2 = sigmoid_vx_obs(vx2)
            for l3 in 1:NV2
                #####
                i3 = TAB_i[l3]
                vx3 = TAB_Vx[i3]
                #
                obs_3 = sigmoid_vx_obs(vx3)
                #####
                m3_obs += TAB_S3[l1, l2, l3] * obs_1 * obs_2 * obs_3
                #####
            end
        end
    end
    #####
    factor = DV^(6) * MASS_PART^2
    #
    m3_obs *= factor
    #####
    return m3_obs
end

##################################################
# Bounded anisotropy observable: O = tanh(vx^2 - vy^2)
##################################################
function tanh_aniso_obs(vx, vy)
    return tanh(vx * vx - vy * vy)
end

function compute_avg_tanh_aniso()
    Obs_val = 0.0
    #####
    for l in 1:NV2
        #####
        i, j = TAB_i[l], TAB_j[l]
        #
        vx = TAB_Vx[i]
        vy = TAB_Vy[j]
        f_val = TAB_F[l]
        #
        Obs_val += f_val * tanh_aniso_obs(vx, vy)
        ####
    end
    Obs_val *= DV^2
    return Obs_val
end

function compute_var_tanh_aniso()
    #####
    Obs_var = 0.0
    #####
    @inbounds for l1 in 1:NV2
        #####
        i1, j1 = TAB_i[l1], TAB_j[l1]
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]
        #
        Obs_1 = tanh_aniso_obs(vx1, vy1)
        #####
        for l2 in 1:NV2
            i2, j2 = TAB_i[l2], TAB_j[l2]
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            #
            Obs_2 = tanh_aniso_obs(vx2, vy2)
            #####
            s2_val = TAB_S2[l1, l2]
            Obs_var += Obs_1 * s2_val * Obs_2
        end
    end
    Obs_var *= DV^(4) * MASS_PART
    #####
    return Obs_var
end

function compute_M3_tanh_aniso()
    #####
    m3_obs = 0.0
    #####
    @inbounds for l1 in 1:NV2
        #####
        i1, j1 = TAB_i[l1], TAB_j[l1]
        #####
        vx1 = TAB_Vx[i1]
        vy1 = TAB_Vy[j1]
        #
        obs_1 = tanh_aniso_obs(vx1, vy1)
        #####
        for l2 in 1:NV2
            #####
            i2, j2 = TAB_i[l2], TAB_j[l2]
            vx2 = TAB_Vx[i2]
            vy2 = TAB_Vy[j2]
            #
            obs_2 = tanh_aniso_obs(vx2, vy2)
            for l3 in 1:NV2
                #####
                i3, j3 = TAB_i[l3], TAB_j[l3]
                vx3 = TAB_Vx[i3]
                vy3 = TAB_Vy[j3]
                #
                obs_3 = tanh_aniso_obs(vx3, vy3)
                #####
                m3_obs += TAB_S3[l1, l2, l3] * obs_1 * obs_2 * obs_3
                #####
            end
        end
    end
    #####
    factor = DV^(6) * MASS_PART^2
    #
    m3_obs *= factor
    #####
    return m3_obs
end

##################################################
# Fused computation of all observables
##################################################

# Each column of TAB_OBSERVABLE_WEIGHTS corresponds to one velocity cell.
# Keeping the observable index first makes the inner loop contiguous in memory.
const OBS_MASS         = 1
const OBS_ENERGY       = 2
const OBS_PX           = 3
const OBS_PY           = 4
const OBS_PTOT         = 5
const OBS_TX           = 6
const OBS_TY           = 7
const OBS_SXXX         = 8
const OBS_SXXY         = 9
const OBS_SXYY         = 10
const OBS_SYYY         = 11
const OBS_NON_POLY     = 12
const OBS_QXXXX        = 13
const OBS_QXXXY        = 14
const OBS_QXXYY        = 15
const OBS_QXYYY        = 16
const OBS_QYYYY        = 17
const OBS_POLYSINCOS   = 18
const OBS_GAUSSIAN     = 19
const OBS_COS_PLUS     = 20
const OBS_COS_MINUS    = 21
const OBS_VSIN_VSIN    = 22
const OBS_TANH_VX      = 23
const OBS_SIGMOID_VX   = 24
const OBS_TANH_ANISO   = 25
const N_OBSERVABLES    = 25

const TAB_OBSERVABLE_WEIGHTS = Matrix{Float64}(undef, N_OBSERVABLES, NV2)

# Padding gives each thread-local accumulator its own cache lines.
const OBS_ACCUMULATOR_STRIDE = 32
const TAB_AVG_OBSERVABLE_ACC = zeros(Float64, N_OBSERVABLES)
const TAB_VAR_OBSERVABLE_ACC = zeros(Float64, OBS_ACCUMULATOR_STRIDE, Threads.nthreads()+1)
const TAB_M3_OBSERVABLE_ACC  = zeros(Float64, OBS_ACCUMULATOR_STRIDE, Threads.nthreads()+1)
const TAB_REDUCED_OBSERVABLE_ACC = zeros(Float64, N_OBSERVABLES)

##################################################
# Precompute all one-particle observables on the velocity grid.
# The elementary observable functions above are deliberately kept as the
# unique definitions of the non-polynomial observables.
##################################################
function initialize_observable_weights!()
    @inbounds for l in 1:NV2
        i, j = TAB_i[l], TAB_j[l]
        vx = TAB_Vx[i]
        vy = TAB_Vy[j]

        vx2 = vx * vx
        vy2 = vy * vy

        TAB_OBSERVABLE_WEIGHTS[OBS_MASS, l]       = 1.0
        TAB_OBSERVABLE_WEIGHTS[OBS_ENERGY, l]     = vx2 + vy2
        TAB_OBSERVABLE_WEIGHTS[OBS_PX, l]         = vx
        TAB_OBSERVABLE_WEIGHTS[OBS_PY, l]         = vy
        TAB_OBSERVABLE_WEIGHTS[OBS_PTOT, l]       = sqrt(vx2 + vy2)
        TAB_OBSERVABLE_WEIGHTS[OBS_TX, l]         = vx2
        TAB_OBSERVABLE_WEIGHTS[OBS_TY, l]         = vy2
        TAB_OBSERVABLE_WEIGHTS[OBS_SXXX, l]       = vx2 * vx
        TAB_OBSERVABLE_WEIGHTS[OBS_SXXY, l]       = vx2 * vy
        TAB_OBSERVABLE_WEIGHTS[OBS_SXYY, l]       = vx * vy2
        TAB_OBSERVABLE_WEIGHTS[OBS_SYYY, l]       = vy2 * vy
        TAB_OBSERVABLE_WEIGHTS[OBS_NON_POLY, l]   = observable_non_polyn(vx, vy)
        TAB_OBSERVABLE_WEIGHTS[OBS_QXXXX, l]      = vx2 * vx2
        TAB_OBSERVABLE_WEIGHTS[OBS_QXXXY, l]      = vx2 * vx * vy
        TAB_OBSERVABLE_WEIGHTS[OBS_QXXYY, l]      = vx2 * vy2
        TAB_OBSERVABLE_WEIGHTS[OBS_QXYYY, l]      = vx * vy2 * vy
        TAB_OBSERVABLE_WEIGHTS[OBS_QYYYY, l]      = vy2 * vy2
        TAB_OBSERVABLE_WEIGHTS[OBS_POLYSINCOS, l] = poly_sin_cos_obs(vx, vy)
        TAB_OBSERVABLE_WEIGHTS[OBS_GAUSSIAN, l]   = gaussian_obs(vx, vy)
        TAB_OBSERVABLE_WEIGHTS[OBS_COS_PLUS, l]   = cos_plus_obs(vx, vy)
        TAB_OBSERVABLE_WEIGHTS[OBS_COS_MINUS, l]  = cos_minus_obs(vx, vy)
        TAB_OBSERVABLE_WEIGHTS[OBS_VSIN_VSIN, l]  = vsin_vsin_obs(vx, vy)
        TAB_OBSERVABLE_WEIGHTS[OBS_TANH_VX, l]    = tanh_vx_obs(vx)
        TAB_OBSERVABLE_WEIGHTS[OBS_SIGMOID_VX, l] = sigmoid_vx_obs(vx)
        TAB_OBSERVABLE_WEIGHTS[OBS_TANH_ANISO, l] = tanh_aniso_obs(vx, vy)
    end
end

initialize_observable_weights!()

##################################################
# Give all three fused computations the same named result layout.
##################################################
@inline function observable_results(accumulator, factor)
    return (
        mass         = accumulator[OBS_MASS]       * factor,
        energy       = accumulator[OBS_ENERGY]     * factor,
        px           = accumulator[OBS_PX]         * factor,
        py           = accumulator[OBS_PY]         * factor,
        ptot         = accumulator[OBS_PTOT]       * factor,
        tx           = accumulator[OBS_TX]         * factor,
        ty           = accumulator[OBS_TY]         * factor,
        sxxx         = accumulator[OBS_SXXX]       * factor,
        sxxy         = accumulator[OBS_SXXY]       * factor,
        sxyy         = accumulator[OBS_SXYY]       * factor,
        syyy         = accumulator[OBS_SYYY]       * factor,
        non_poly     = accumulator[OBS_NON_POLY]   * factor,
        qxxxx        = accumulator[OBS_QXXXX]      * factor,
        qxxxy        = accumulator[OBS_QXXXY]      * factor,
        qxxyy        = accumulator[OBS_QXXYY]      * factor,
        qxyyy        = accumulator[OBS_QXYYY]      * factor,
        qyyyy        = accumulator[OBS_QYYYY]      * factor,
        poly_sin_cos = accumulator[OBS_POLYSINCOS] * factor,
        gaussian     = accumulator[OBS_GAUSSIAN]   * factor,
        cos_plus     = accumulator[OBS_COS_PLUS]   * factor,
        cos_minus    = accumulator[OBS_COS_MINUS]  * factor,
        vsin_vsin    = accumulator[OBS_VSIN_VSIN]  * factor,
        tanh_vx      = accumulator[OBS_TANH_VX]    * factor,
        sigmoid_vx   = accumulator[OBS_SIGMOID_VX] * factor,
        tanh_aniso   = accumulator[OBS_TANH_ANISO] * factor,
    )
end

##################################################
# All averages: one traversal of F.
##################################################
function compute_all_averages()
    fill!(TAB_AVG_OBSERVABLE_ACC, 0.0)

    @inbounds for l in 1:NV2
        f_val = TAB_F[l]
        @simd for iobs in 1:N_OBSERVABLES
            TAB_AVG_OBSERVABLE_ACC[iobs] = muladd(
                f_val,
                TAB_OBSERVABLE_WEIGHTS[iobs, l],
                TAB_AVG_OBSERVABLE_ACC[iobs],
            )
        end
    end

    return observable_results(TAB_AVG_OBSERVABLE_ACC, DV^2)
end

##################################################
# All variances: one cache-friendly, threaded traversal of S2.
##################################################
function compute_all_variances()
    fill!(TAB_VAR_OBSERVABLE_ACC, 0.0)

    Threads.@threads :static for l2 in 1:NV2
        tid = Threads.threadid()
        @inbounds for l1 in 1:NV2
            s2_val = TAB_S2[l1, l2]
            @simd for iobs in 1:N_OBSERVABLES
                weight = TAB_OBSERVABLE_WEIGHTS[iobs, l1] * TAB_OBSERVABLE_WEIGHTS[iobs, l2]
                TAB_VAR_OBSERVABLE_ACC[iobs, tid] = muladd(
                    s2_val,
                    weight,
                    TAB_VAR_OBSERVABLE_ACC[iobs, tid],
                )
            end
        end
    end

    fill!(TAB_REDUCED_OBSERVABLE_ACC, 0.0)
    @inbounds for tid in 1:Threads.nthreads()
        @simd for iobs in 1:N_OBSERVABLES
            TAB_REDUCED_OBSERVABLE_ACC[iobs] += TAB_VAR_OBSERVABLE_ACC[iobs, tid]
        end
    end

    return observable_results(TAB_REDUCED_OBSERVABLE_ACC, DV^4 * MASS_PART)
end

##################################################
# All third moments: one cache-friendly, threaded traversal of S3.
##################################################
function compute_all_M3()
    fill!(TAB_M3_OBSERVABLE_ACC, 0.0)

    Threads.@threads :static for l3 in 1:NV2
        tid = Threads.threadid()
        @inbounds for l2 in 1:NV2
            for l1 in 1:NV2
                s3_val = TAB_S3[l1, l2, l3]
                @simd for iobs in 1:N_OBSERVABLES
                    weight = TAB_OBSERVABLE_WEIGHTS[iobs, l1] *
                             TAB_OBSERVABLE_WEIGHTS[iobs, l2] *
                             TAB_OBSERVABLE_WEIGHTS[iobs, l3]
                    TAB_M3_OBSERVABLE_ACC[iobs, tid] = muladd(
                        s3_val,
                        weight,
                        TAB_M3_OBSERVABLE_ACC[iobs, tid],
                    )
                end
            end
        end
    end

    fill!(TAB_REDUCED_OBSERVABLE_ACC, 0.0)
    @inbounds for tid in 1:Threads.nthreads()
        @simd for iobs in 1:N_OBSERVABLES
            TAB_REDUCED_OBSERVABLE_ACC[iobs] += TAB_M3_OBSERVABLE_ACC[iobs, tid]
        end
    end

    return observable_results(TAB_REDUCED_OBSERVABLE_ACC, DV^6 * MASS_PART^2)
end

##################################################
# Store fused observable results in the output arrays.
##################################################
function store_averages!(index)
    avg = compute_all_averages()

    TAB_E[index]  = avg.energy
    TAB_M[index]  = avg.mass
    TAB_Px[index] = avg.px
    TAB_Py[index] = avg.py
    TAB_Tx[index] = avg.tx
    TAB_Ty[index] = avg.ty

    TAB_SXXX[index] = avg.sxxx
    TAB_SXXY[index] = avg.sxxy
    TAB_SXYY[index] = avg.sxyy
    TAB_SYYY[index] = avg.syyy
    TAB_NON_POLY_OBS[index] = avg.non_poly

    TAB_QXXXX[index] = avg.qxxxx
    TAB_QXXXY[index] = avg.qxxxy
    TAB_QXXYY[index] = avg.qxxyy
    TAB_QXYYY[index] = avg.qxyyy
    TAB_QYYYY[index] = avg.qyyyy

    TAB_POLYSINCOS[index] = avg.poly_sin_cos
    TAB_GAUSSIAN[index]   = avg.gaussian
    TAB_COS_PLUS[index]   = avg.cos_plus
    TAB_COS_MINUS[index]  = avg.cos_minus
    TAB_VSIN_VSIN[index]  = avg.vsin_vsin
    TAB_TANH_VX[index]    = avg.tanh_vx
    TAB_SIGMOID_VX[index] = avg.sigmoid_vx
    TAB_TANH_ANISO[index] = avg.tanh_aniso
end

function store_variances!(index)
    var = compute_all_variances()

    TAB_VAR_M[index]    = var.mass
    TAB_VAR_E[index]    = var.energy
    TAB_VAR_Px[index]   = var.px
    TAB_VAR_Py[index]   = var.py
    TAB_VAR_Ptot[index] = var.ptot
    TAB_VAR_Tx[index]   = var.tx
    TAB_VAR_Ty[index]   = var.ty

    TAB_VAR_SXXX[index] = var.sxxx
    TAB_VAR_SXXY[index] = var.sxxy
    TAB_VAR_SXYY[index] = var.sxyy
    TAB_VAR_SYYY[index] = var.syyy
    TAB_VAR_NON_POLY_OBS[index] = var.non_poly

    TAB_VAR_QXXXX[index] = var.qxxxx
    TAB_VAR_QXXXY[index] = var.qxxxy
    TAB_VAR_QXXYY[index] = var.qxxyy
    TAB_VAR_QXYYY[index] = var.qxyyy
    TAB_VAR_QYYYY[index] = var.qyyyy

    TAB_VAR_POLYSINCOS[index] = var.poly_sin_cos
    TAB_VAR_GAUSSIAN[index]   = var.gaussian
    TAB_VAR_COS_PLUS[index]   = var.cos_plus
    TAB_VAR_COS_MINUS[index]  = var.cos_minus
    TAB_VAR_VSIN_VSIN[index]  = var.vsin_vsin
    TAB_VAR_TANH_VX[index]    = var.tanh_vx
    TAB_VAR_SIGMOID_VX[index] = var.sigmoid_vx
    TAB_VAR_TANH_ANISO[index] = var.tanh_aniso
end

function store_M3!(index)
    m3 = compute_all_M3()

    TAB_M3_M[index]    = m3.mass
    TAB_M3_E[index]    = m3.energy
    TAB_M3_Px[index]   = m3.px
    TAB_M3_Py[index]   = m3.py
    TAB_M3_Ptot[index] = m3.ptot
    TAB_M3_Tx[index]   = m3.tx
    TAB_M3_Ty[index]   = m3.ty

    TAB_M3_SXXX[index] = m3.sxxx
    TAB_M3_SXXY[index] = m3.sxxy
    TAB_M3_SXYY[index] = m3.sxyy
    TAB_M3_SYYY[index] = m3.syyy
    TAB_M3_NON_POLY_OBS[index] = m3.non_poly

    TAB_M3_QXXXX[index] = m3.qxxxx
    TAB_M3_QXXXY[index] = m3.qxxxy
    TAB_M3_QXXYY[index] = m3.qxxyy
    TAB_M3_QXYYY[index] = m3.qxyyy
    TAB_M3_QYYYY[index] = m3.qyyyy

    TAB_M3_POLYSINCOS[index] = m3.poly_sin_cos
    TAB_M3_GAUSSIAN[index]   = m3.gaussian
    TAB_M3_COS_PLUS[index]   = m3.cos_plus
    TAB_M3_COS_MINUS[index]  = m3.cos_minus
    TAB_M3_VSIN_VSIN[index]  = m3.vsin_vsin
    TAB_M3_TANH_VX[index]    = m3.tanh_vx
    TAB_M3_SIGMOID_VX[index] = m3.sigmoid_vx
    TAB_M3_TANH_ANISO[index] = m3.tanh_aniso
end

function store_observables!(index)
    store_averages!(index)
    if HIERARCHY == "S2"
        store_variances!(index)
    elseif HIERARCHY == "S3"
        store_M3!(index)
    end
end

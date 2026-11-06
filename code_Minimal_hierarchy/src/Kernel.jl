##################################################
# Compute all values of B(v, vp) in the NV x NV grid
# TAB_B_{apha,beta}[i,j][ip,jp] = B(vx[i] - vx'[ip], vy[j] - vy'[jp])
# @IMPROVE -- Should put parallelisation
##################################################
function TAB_B!()
    Threads.@threads for lp in 1:NV2
        #####
        ip, jp = TAB_i[lp], TAB_j[lp]
        vxp = TAB_Vx[ip]
        vyp = TAB_Vy[jp]
        #####
        for l in 1:NV2
            i, j = TAB_i[l], TAB_j[l]
            vy = TAB_Vy[j]
            vx = TAB_Vx[i]
            #####
            ux = vx - vxp
            uy = vy - vyp
            #####
            Bxx, Bxy, Byx, Byy = get_kernel_B(ux,uy)
            #####
            TAB_Bxx[l,lp] = Bxx
            TAB_Bxy[l,lp] = Bxy
            TAB_Byy[l,lp] = Byy
            #####
        end
    end
end

##################################################
# Initialising the collision kernel, once and for all
TAB_B!()
##################################################
const TAB_TAB_Bxx = [ [TAB_Bxx[l, lp] for lp=1:NV2] for l=1:NV2]
const TAB_TAB_Bxy = [ [TAB_Bxy[l, lp] for lp=1:NV2] for l=1:NV2]
const TAB_TAB_Byy = [ [TAB_Byy[l, lp] for lp=1:NV2] for l=1:NV2]

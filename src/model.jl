using Printf
using Dates
using NCDatasets
using SparseArrays, LinearAlgebra


#make sum drop dimension automatically and set subnormals to zero to avoid slowdown
sumd(x,dims) = dropdims(sum(x,dims=dims),dims=dims)
set_zero_subnormals(true)


#TODO add a sinking (dsink) to remove some detritus or it wont equilibriate 
function run_NPZBDV(prms, lysis, pulse=0)

    start_time = now()
    trec = prms.nt ÷ prms.nrec 

    # Generate empty arrays
    nrec1 = Int(prms.nrec + 1)
    track_n = Array{Float64,2}(undef, prms.nn, nrec1) 
    track_c = Array{Float64,2}(undef, prms.nc, nrec1) 
    track_p = Array{Float64,2}(undef, prms.np, nrec1) 
    track_z = Array{Float64,2}(undef, prms.nz, nrec1) 
    track_b = Array{Float64,2}(undef, prms.nb, nrec1) 
    track_dn = Array{Float64,2}(undef, prms.ndn, nrec1)
    track_dc = Array{Float64,2}(undef, prms.ndc, nrec1)
    track_v = Array{Float64,2}(undef, prms.nv, nrec1)
    track_time = Array{Float64,1}(undef, nrec1)  

    track_n[:,1] .= prms.nIC
    track_c[:,1] .= prms.cIC
    track_p[:,1] .= prms.pIC
    track_z[:,1] .= prms.zIC
    track_b[:,1] .= prms.bIC
    track_dn[:,1] .= prms.dnIC
    track_dc[:,1] .= prms.dcIC
    track_v[:,1] .= prms.vIC

    track_time[1] = 0

    #--------------------------------------
    #Initial conditions at time t = 0
    ptemp = copy(prms.pIC) 
    btemp = copy(prms.bIC)
    ztemp = copy(prms.zIC) 
    ntemp = copy(prms.nIC) 
    ctemp = copy(prms.cIC) 
    dntemp = copy(prms.dnIC) 
    dctemp = copy(prms.dcIC)
    vtemp = copy(prms.vIC) 


    for t = 1:prms.nt 

        # Solver - Runge-Kutta 4th order 
        ntemp, ctemp, ptemp, ztemp, btemp, dntemp, dctemp, vtemp = rk4(ntemp, ctemp, ptemp, ztemp, btemp, dntemp, dctemp, vtemp, prms, t, lysis)

        if mod(t, trec)==0

            track_n,track_c,track_p,track_z,track_b,track_dn,track_dc,track_v,track_time = update_tracking_arrs(track_n, track_c, track_p, track_z, track_b, track_dn, track_dc, track_v, track_time, 
                                                                                ntemp, ctemp, ptemp, ztemp, btemp, dntemp, dctemp, vtemp, t, trec, prms)
            
            if prms.carbon == 1
                tot_n = sum(ntemp) + sum(dntemp) + ((sum(ptemp) + sum(btemp) + sum(ztemp) + sum(vtemp)) * (1/prms.CNr))
                tot_c = sum(ctemp) + sum(dctemp) + ((sum(ptemp) + sum(btemp) + sum(ztemp) + sum(vtemp)) * (1-(1/prms.CNr)))
                println("Total N: ", tot_n)
                println("Total C: ", tot_c)
            else
                tot_n = sum(ntemp) + sum(dntemp) + sum(ptemp) + sum(btemp) + sum(ztemp) + sum(vtemp)
                println("Total N: ", tot_n)
            end


        end 

        # if pulse != 0
        #     if prms.pulse == 1
        #         # 1 = winter, pulse 5.0 every 10 days
        #         if t % 1000 == 0
        #             pulse = nutrient_pulse(prms.pulse)
        #             dtemp .+= pulse     
        #             ntemp .+= pulse            
        #         end
        #     else 
        #         # 2 = summer, pulse 2.0 every 30 days
        #         if t % 3000 == 0
        #             pulse = nutrient_pulse(prms.pulse)
        #             dtemp .+= pulse     
        #             ntemp .+= pulse        
        #         end
        #     end
        # end

        if t == nt
            
            end_time = now() 
            savetoNC(fsaven, track_p, track_b, track_z, track_n, track_c, track_dn, track_dc, track_v, track_time, start_time, end_time, prms, pulse)

        end
    end 


    return ntemp, ctemp, ptemp, ztemp, btemp, dntemp, dctemp, vtemp, track_time, fsaven

end 


function model_functions(N, C, P, Z, B, Dn, Dc, V, prms, t, lysis)

    dNdt = zeros(Float64, prms.nn)
    dCdt = zeros(Float64, prms.nc)
    dPdt = zeros(Float64, prms.np)
    dZdt = zeros(Float64, prms.nz)
    dBdt = zeros(Float64, prms.nb)
    dDndt = zeros(Float64, prms.ndn)
    dDcdt = zeros(Float64, prms.ndc)
    dVdt = zeros(Float64, prms.nv)

    don_gain_mort = zeros(Float64, 1)
    don_gain_vly = zeros(Float64, 1)
    don_gain_vde = zeros(Float64, 1)

    doc_gain_mort = zeros(Float64, 1)
    doc_gain_vly = zeros(Float64, 1)
    doc_gain_vde = zeros(Float64, 1)

    # DIN and DIC supply (turn off rsource and rsink to check if conserving)
    if prms.carbon == 1
        dNdt .+= prms.rsource * (1/prms.CNr)
        dCdt .+= prms.rsource * (1-(1/prms.CNr))
    else
        dNdt .+= prms.rsource
    end

    # phyto uptake 
    dPdt, dNdt, dCdt = phyto_uptake(prms, N, C, P, dNdt, dCdt, dPdt, t)

    # bacteria uptake
    dDndt, dDcdt, dBdt, dNdt, dCdt = bacteria_uptake(prms, B, Dn, Dc, dDndt, dDcdt, dBdt, dNdt, dCdt, t)

    #phytoplankton mortality
    dPdt, don_gain_mort, doc_gain_mort = phyto_mortality(prms, P, dPdt, don_gain_mort, doc_gain_mort, t)

    #bacterial mortality
    dBdt, don_gain_mort, doc_gain_mort = bacterial_mortality(prms, B, dBdt, don_gain_mort, doc_gain_mort, lysis, t)
    
    if prms.graze == 1
        # zooplank grazing
        dZdt, dNdt, dCdt, dPdt, dBdt = grazing(prms, P, B, Z, dZdt, dNdt, dCdt, dPdt, dBdt, t)
        #zooplankton mortality 
        dZdt, don_gain_mort, doc_gain_mort = zoo_mortality(prms, Z, dZdt, don_gain_mort, doc_gain_mort, t)
    else
    end

    if prms.lysis == 1
        # viral lysis (B only for now)
        dVdt, dBdt, don_gain_vly, doc_gain_vly = viral_lysis(prms, B, V, dVdt, dBdt, don_gain_vly, doc_gain_vly, t)
        #viral decay
        dVdt, don_gain_vde, doc_gain_vde = viral_decay(prms, V, dVdt, don_gain_vde, doc_gain_vde, t)
    else
    end

    #sinking rate and split accumulated OM into nd pools
    dDndt, dDcdt = total_change_in_d(prms, Dn, Dc, dDndt, dDcdt, don_gain_mort, doc_gain_mort, don_gain_vly, doc_gain_vly, don_gain_vde, doc_gain_vde, t)

    return dNdt, dCdt, dPdt, dZdt, dBdt, dDndt, dDcdt, dVdt

end 


function phyto_uptake(prms, N, C, P, dNdt, dCdt, dPdt, t)
    # NOTE phyto consume 6 C for every 1 N 
    # NOTE this could be set to scale according to environmental C:N ratio

    II, JJ = get_nonzero_axes(prms.CMp)

    for j = axes(II, 1)

        if prms.carbon == 1
            # uptake rate set by limiting nutrient (C or N)
            mu = prms.vmax_ij[II[j],JJ[j]] .* min.(N ./ (N .+ prms.Kp_ij[II[j],JJ[j]]), C ./ (C .+ prms.Kp_ij[II[j],JJ[j]]))
            uptake = P[JJ[j],:] .* mu
            dPdt[JJ[j],:] += uptake
            dNdt -= uptake * (1/prms.CNr)
            dCdt -= uptake * (1-(1/prms.CNr))
        else
            mu = prms.vmax_ij[II[j],JJ[j]] .* (N ./ (N .+ prms.Kp_ij[II[j],JJ[j]]))
            uptake = P[JJ[j],:] .* mu
            dPdt[JJ[j],:] += uptake
            dNdt -= uptake
        end
    end

    return dPdt, dNdt, dCdt

end 


function bacteria_uptake(prms, B, Dn, Dc, dDndt, dDcdt, dBdt, dNdt, dCdt, t=0)
    #NOTE for every 1 N excreted, CNr * C is respired
    # uptake rate controlled by limiting substrate (DOC or DON)

    II, JJ = get_nonzero_axes(prms.CM)

    for j = axes(II, 1)
        # DON and DOC uptake and remineralization (1 - yield)
        if prms.carbon == 1
            mu = prms.umax_ij[II[j],JJ[j]] .* min.( Dn[II[j],:]./(Dn[II[j],:].+prms.Km_ij[II[j],JJ[j]]), Dc[II[j],:]./(Dc[II[j],:].+prms.Km_ij[II[j],JJ[j]]) )
            uptake = B[JJ[j],:] .* mu
            dBdt[JJ[j],:] += (uptake * prms.y_ij[II[j],JJ[j]]) 
            dDndt[II[j],:] -= (uptake * (1/prms.CNr))
            dDcdt[II[j],:] -= (uptake * (1-(1/prms.CNr)))

            # proportion remineralized
            dNdt += (uptake * (1-prms.y_ij[II[j],JJ[j]])) * (1/prms.CNr)
            dCdt += (uptake * (1-prms.y_ij[II[j],JJ[j]])) * (1-(1/prms.CNr))
        else 
            mu = prms.umax_ij[II[j],JJ[j]] .* Dn[II[j],:]./(Dn[II[j],:].+prms.Km_ij[II[j],JJ[j]])
            uptake = B[JJ[j],:] .* mu
            dBdt[JJ[j],:] += (uptake * prms.y_ij[II[j],JJ[j]]) 
            dDndt[II[j],:] -= uptake
            dNdt += uptake * (1-prms.y_ij[II[j],JJ[j]]) 
        end

    end

    return dDndt, dDcdt, dBdt, dNdt, dCdt

end


function grazing(prms, P, B, Z, dZdt, dNdt, dCdt, dPdt, dBdt, t=0)

    GrM = copy(prms.GrM)
    for k = 1:prms.nz
        if sum(GrM[k, 1:prms.np]) > 0 
            dZdt, dNdt, dCdt, dPdt = phyto_grazing(prms, GrM, P, Z, dZdt, dNdt, dCdt, dPdt, k, t)
        end
        if sum(GrM[k, prms.np+1:end]) > 0 
            dZdt, dNdt, dCdt, dBdt = bacteria_grazing(prms, GrM, B, Z, dZdt, dNdt, dCdt, dBdt, k, t)
        end
    end

    return dZdt, dNdt, dCdt, dPdt, dBdt

end

        function phyto_grazing(prms, GrM, P, Z, dZdt, dNdt, dCdt, dPdt, k, t)

            prey = sum(GrM[k,1:prms.np] .* P)
            gz = prms.g_max[k] * prey / (prey + prms.K_g[k])
            dZdt[k,:] += prms.γ[k] * gz .* Z[k,:]
            dPdt -= (gz .* Z[k,:] .* GrM[k,1:prms.np] .* P ./ prey)

            if prms.carbon == 1
                dNdt += ((1-prms.γ[k]) .* gz .* Z[k,:]) * (1/prms.CNr)
                dCdt += ((1-prms.γ[k]) .* gz .* Z[k,:]) * (1-(1/prms.CNr))
            else
                dNdt += ((1-prms.γ[k]) .* gz .* Z[k,:]) 
            end
            
            return dZdt, dNdt, dCdt, dPdt

        end

        function bacteria_grazing(prms, GrM, B, Z, dZdt, dNdt, dCdt, dBdt, k, t)

            prey = sum(GrM[k,prms.np+1:end] .* B)
            gz = prms.g_max[k] * prey / (prey + prms.K_g[k])
            dZdt[k,:] += prms.γ[k] .* gz .* Z[k,:]
            dBdt -=  (gz .* Z[k,:] .* GrM[k,prms.np+1:end] .* B ./ prey)

        if prms.carbon == 1
            dNdt += ((1 - prms.γ[k]) .* gz .* Z[k,:]) * (1/prms.CNr)
            dCdt += ((1 - prms.γ[k]) .* gz .* Z[k,:]) * (1-(1/prms.CNr))
        else
            dNdt += ((1 - prms.γ[k]) .* gz .* Z[k,:]) 
        end

            return dZdt, dNdt, dCdt, dBdt

        end


function viral_lysis(prms, B, V, dVdt, dBdt, don_gain_vly, doc_gain_vly, t)

    Qb = 2e-12             # mmol N/cell  (estimate from Weitz et al range of 0.5 to 4e-12 mmol N / cell)
    Qv = 1e-14             # mmol N/cell (Weitz et al 2015 range is 1.4 to 14e-15 mmol N / cell)
    lysis_rate = 8e-15     # m3/virus/day
    burst_size = 25        # 1 infected B cell = 25 viruses

    II, JJ = get_nonzero_axes(prms.VM)
    for j = axes(II, 1)
        # lysis_Bi = VM[II[j], JJ[j]] .* lysis_rate .* (V[II[j],:] * 1/Qv)  .* (B[JJ[j],:] * 1/Qb)    # cells/m3  (original weitz version)
        lysis_Bi = VM[II[j], JJ[j]] .* lysis_rate .* (V[II[j],:] * 1/Qv)  .* ((B[JJ[j],:] .* B[JJ[j],:]) * 1/Qb)     # cells/m3 quadratic version
        bacteria_loss = lysis_Bi .* Qb                                                              # cells/m3 * mmol N/cell = mmol N/m3
        viral_growth = lysis_Bi .* burst_size .* Qv                                                 # cells/m3 * 25 * mmol N/cell = mmol N/m3
        dVdt[II[j],:] += viral_growth
        dBdt[JJ[j],:] -= bacteria_loss


        released_to_dom = bacteria_loss - viral_growth
        if prms.carbon == 1
            don_gain_vly += released_to_dom * (1/prms.CNr)
            doc_gain_vly += released_to_dom * (1-(1/prms.CNr))
        else
            don_gain_vly += released_to_dom
        end

    end

    return dVdt, dBdt, don_gain_vly, doc_gain_vly

end


function phyto_mortality(prms, P, dPdt, don_gain_mort, doc_gain_mort, t=0)

    pmort = (prms.m_lp .+ prms.m_qp .* P) .* P
    dPdt -= pmort

    if prms.carbon == 1
        don_gain_mort .+= sum(pmort) * (1/prms.CNr)
        doc_gain_mort .+= sum(pmort) * (1-(1/prms.CNr))
    else
        don_gain_mort .+= sum(pmort)
    end

    return dPdt, don_gain_mort, doc_gain_mort

end


function bacterial_mortality(prms, B, dBdt, don_gain_mort, doc_gain_mort, lysis, t=0)

    if prms.lysis == 1
        bmort = prms.m_lb .* B
    else
        bmort = (prms.m_lb .+ prms.m_qb .* B) .* B
    end

    dBdt -= bmort

    if prms.carbon == 1
        don_gain_mort .+= sum(bmort) * (1/prms.CNr)
        doc_gain_mort .+= sum(bmort) * (1-(1/prms.CNr))
    else
        don_gain_mort .+= sum(bmort)
    end

    return dBdt, don_gain_mort, doc_gain_mort

end


function zoo_mortality(prms, Z, dZdt, don_gain_mort, doc_gain_mort, t=0)

    zmort = (prms.m_lz .+ prms.m_qz .* Z) .* Z
    dZdt -= zmort

    if prms.carbon == 1
        don_gain_mort .+= sum(zmort) * (1/prms.CNr)
        doc_gain_mort .+= sum(zmort) * (1-(1/prms.CNr))
    else 
        don_gain_mort .+= sum(zmort)
    end

    return dZdt, don_gain_mort, doc_gain_mort

end


function viral_decay(prms, V, dVdt, don_gain_vde, doc_gain_vde, t=0)

    decay = prms.vde .* V
    dVdt -= decay

    if prms.carbon == 1
        don_gain_vde .+= sum(decay) * (1/prms.CNr)
        doc_gain_vde .+= sum(decay) * (1-(1/prms.CNr))
    else
        don_gain_vde .+= sum(decay)
    end

    return dVdt, don_gain_vde, doc_gain_vde
end


function total_change_in_d(prms, Dn, Dc, dDndt, dDcdt, don_gain_mort, doc_gain_mort, don_gain_vly, doc_gain_vly, don_gain_vde, doc_gain_vde, t=0)

    dDndt += don_gain_mort .* prms.om_dist_mort
    dDndt += don_gain_vly .* prms.om_dist_lys
    dDndt += don_gain_vde .* prms.om_dist_vde
    dDndt -= Dn .* prms.rsink

    if prms.carbon == 1
        dDcdt += doc_gain_mort .* prms.om_dist_mort
        dDcdt += doc_gain_vly .* prms.om_dist_lys
        dDcdt += doc_gain_vde .* prms.om_dist_vde
        dDcdt -= Dc .* prms.rsink 
    end

    return dDndt, dDcdt

end


function get_nonzero_axes(M)

    Cs = sparse(M)
    (II, JJ, _) = findnz(Cs) 
    
    return II, JJ

end 


function nutrient_pulse(pulse)

    if pulse == 1
        pulse_size = 5.0
    else
        pulse_size = 2.0
    end

    return pulse_size

end

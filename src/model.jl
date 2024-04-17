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
    track_p = Array{Float64,2}(undef, prms.np, nrec1) 
    track_z = Array{Float64,2}(undef, prms.nz, nrec1) 
    track_b = Array{Float64,2}(undef, prms.nb, nrec1) 
    track_d = Array{Float64,2}(undef, prms.nd, nrec1)
    track_v = Array{Float64,2}(undef, prms.nv, nrec1)
    track_time = Array{Float64,1}(undef, nrec1)  

    track_n[:,1] .= prms.nIC
    track_p[:,1] .= prms.pIC
    track_z[:,1] .= prms.zIC
    track_b[:,1] .= prms.bIC
    track_d[:,1] .= prms.dIC
    track_v[:,1] .= prms.vIC

    track_time[1] = 0

    #--------------------------------------
    #Initial conditions at time t = 0
    ptemp = copy(prms.pIC) 
    btemp = copy(prms.bIC)
    ztemp = copy(prms.zIC) 
    ntemp = copy(prms.nIC) 
    dtemp = copy(prms.dIC) 
    vtemp = copy(prms.vIC) 


    for t = 1:prms.nt 

        # Solver - Runge-Kutta 4th order 
        ntemp, ptemp, ztemp, btemp, dtemp, vtemp = rk4(ntemp, ptemp, ztemp, btemp, dtemp, vtemp, prms, t, lysis)

        if mod(t, trec)==0

            track_n, track_p, track_z, track_b, track_d, track_v, track_time = update_tracking_arrs(track_n, track_p, track_z, track_b, track_d, track_v, track_time, 
                                                                                ntemp, ptemp, ztemp, btemp, dtemp, vtemp, t, trec, prms)
            println("Total N: ", sum(ptemp) + sum(btemp) + sum(ntemp) + sum(dtemp) + sum(ztemp) + sum(vtemp))
        
        end 

        if pulse != 0
            if prms.pulse == 1
                # 1 = winter, pulse 5.0 every 10 days
                if t % 1000 == 0
                    pulse = nutrient_pulse(prms.pulse)
                    dtemp .+= pulse     
                    ntemp .+= pulse            
                end
            else 
                # 2 = summer, pulse 2.0 every 30 days
                if t % 3000 == 0
                    pulse = nutrient_pulse(prms.pulse)
                    dtemp .+= pulse     
                    ntemp .+= pulse        
                end
            end
        end

        #calculate uptake and v for last timepoint
        if t == nt

            # II, JJ = get_nonzero_axes(prms)
            # v = zeros(prms.nd, prms.nb) 
            # uptake = zeros(prms.nd, prms.nb) 
        
            # @inbounds for n = axes(II, 1)
            #     v[II[n],JJ[n]] = prms.vmax_ij[II[n],JJ[n]] * prms.pen[JJ[n]] * dtemp[II[n]]/(dtemp[II[n]] + prms.Km_ij[II[n],JJ[n]]) 
            #     uptake[II[n],JJ[n]] = v[II[n],JJ[n]] * btemp[JJ[n]]
            # end  
            
            end_time = now() 

            savetoNC(fsaven, track_p, track_b, track_z, track_n, track_d, track_v, track_time, start_time, end_time, prms, pulse)

        end
    end 


    return ntemp, ptemp, ztemp, btemp, dtemp, vtemp, track_time, fsaven

end 


function model_functions(N, P, Z, B, D, V, prms, t, lysis)

    dNdt = zeros(Float64, prms.nn)
    dPdt = zeros(Float64, prms.np)
    dZdt = zeros(Float64, prms.nz)
    dBdt = zeros(Float64, prms.nb)
    dDdt = zeros(Float64, prms.nd)
    dVdt = zeros(Float64, prms.nv)

    d_gain_mort = zeros(Float64, 1)
    d_gain_vly = zeros(Float64, 1)
    d_gain_vde = zeros(Float64, 1)

    # OM supply (turn off rsource and rsink to check if N conserving)
    dNdt .+= prms.rsource

    # phyto uptake >>> TODO phyto assimilation efficiency and contribution to dDdt?
    dPdt, dNdt = phyto_uptake(prms, N, P, dNdt, dPdt, t)

    # bacteria uptake
    dDdt, dBdt, dNdt = bacteria_uptake(prms, B, D, dDdt, dBdt, dNdt, t)

    # zooplank grazing
    dZdt, dNdt, dPdt, dBdt = grazing(prms, P, B, Z, dZdt, dNdt, dPdt, dBdt, t)

    #phytoplankton mortality
    dPdt, d_gain_mort = phyto_mortality(prms, P, dPdt, d_gain_mort, t)

    #bacterial mortality
    dBdt, d_gain_mort = bacterial_mortality(prms, B, dBdt, d_gain_mort, lysis, t)
    
    #zooplankton mortality 
    dZdt, d_gain_mort = zoo_mortality(prms, Z, dZdt, d_gain_mort, t)

    if lysis == 1
        # viral lysis (B only for now)
        dVdt, dBdt, d_gain_vly = viral_lysis(prms, B, V, dVdt, dBdt, d_gain_vly, t)
        #viral decay
        dVdt, d_gain_vde = viral_decay(prms, V, dVdt, d_gain_vde, t)
    else
    end

    #sinking rate and split accumulated OM into nd pools
    dDdt = total_change_in_d(prms, D, dDdt, d_gain_mort, d_gain_vly, d_gain_vde, t)

    return dNdt, dPdt, dZdt, dBdt, dDdt, dVdt

end 


function phyto_uptake(prms, N, P, dNdt, dPdt, t)

    II, JJ = get_nonzero_axes(prms.CMp)

    for j = axes(II, 1)
        uptake = P[JJ[j]] .* prms.vmax_ij[II[j],JJ[j]] .* N ./ (N .+ prms.Kp_ij[II[j],JJ[j]])
        dNdt .+= -uptake[1]
        dPdt[JJ[j]] += uptake[1]
    end

    return dPdt, dNdt

end 


function bacteria_uptake(prms, B, D, dDdt, dBdt, dNdt, t=0)

    II, JJ = get_nonzero_axes(prms.CM)

    for j = axes(II, 1)
        uptake = prms.umax_ij[II[j],JJ[j]] .* D[II[j],:] ./ (D[II[j],:] .+ prms.Km_ij[II[j],JJ[j]]) .* B[JJ[j],:]
        dDdt[II[j],:] += -uptake
        dBdt[JJ[j],:] += uptake .* prms.y_ij[II[j],JJ[j]]
        dNdt += uptake * (1 - prms.y_ij[II[j],JJ[j]])
    end

    return dDdt, dBdt, dNdt

end


function grazing(prms, P, B, Z, dZdt, dNdt, dPdt, dBdt, t=0)

    GrM = copy(prms.GrM)
    for k = 1:prms.nz
        if sum(GrM[k, 1:prms.np]) > 0 
            dZdt, dNdt, dPdt = phyto_grazing(prms, GrM, P, Z, dZdt, dNdt, dPdt, k, t)
        end
        if sum(GrM[k, prms.np+1:end]) > 0 
            dZdt, dNdt, dBdt = bacteria_grazing(prms, GrM, B, Z, dZdt, dNdt, dBdt, k, t)
        end
    end

    return dZdt, dNdt, dPdt, dBdt

end

        function phyto_grazing(prms, GrM, P, Z, dZdt, dNdt, dPdt, k, t)

            prey = sum(GrM[k,1:prms.np] .*P)
            gp = prms.g_max[k] * prey / (prey + prms.K_g[k])
            dZdt[k,:] += prms.γ[k] * gp .* Z[k,:]
            dNdt += (1 - prms.γ[k]) .* gp .* Z[k,:]
            dPdt += -gp .* Z[k,:] .* GrM[k,1:prms.np] .* P ./ prey 
            
            return dZdt, dNdt, dPdt

        end

        function bacteria_grazing(prms, GrM, B, Z, dZdt, dNdt, dBdt, k, t)

            prey = sum(GrM[k,prms.np+1:end] .* B)
            gb = prms.g_max[k] * prey / (prey + prms.K_g[k])
            dZdt[k,:] += prms.γ[k] .* gb .* Z[k,:]
            dNdt += (1 - prms.γ[k]) .* gb .* Z[k,:]
            dBdt +=  -gb .* Z[k,:] .* GrM[k,prms.np+1:end] .* B ./ prey

            return dZdt, dNdt, dBdt

        end


function viral_lysis(prms, B, V, dVdt, dBdt, d_gain_vly, t)

    # II, JJ = get_nonzero_axes(prms.VM)

    # for j = axes(II, 1)
    #     lysis_Bi = prms.vly .* VM[II[j], JJ[j]] .* B[JJ[j],:] .* V[II[j],:]
    #     dBdt[JJ[j],:] += -lysis_Bi
    #     dVdt[II[j],:] += lysis_Bi .* prms.vbs 
    #     # dNdt += -(lysis_Bi .* prms.vbs)
    #     d_gain_vly += lysis_Bi
    # end
    #NOTE as a first approximation, 30% of the lysed B goes into V growth, and 70% is returned to D 
    # leaving out burst size for now as it was just killing all B fast
    # lysis = prms.vly .* B .* V .* prms.vbs
    lysis = prms.vly .* B .* V 
    v_growth = lysis .* 0.4
    d_gain_vly += [sum(lysis) * 0.6]
    dVdt += [sum(v_growth)]
    dBdt += -lysis

    return dVdt, dBdt, d_gain_vly

end


function phyto_mortality(prms, P, dPdt, d_gain_mort=0, t=0)

    pmort = (prms.m_lp .+ prms.m_qp .* P) .* P
    dPdt += -pmort
    d_gain_mort .+= sum(pmort)

    return dPdt, d_gain_mort

end


function bacterial_mortality(prms, B, dBdt, d_gain_mort, lysis, t=0)

    if lysis == 0
        bmort = (prms.m_lb .+ prms.m_qb .* B) .* B
    else
        bmort = prms.m_lb .* B
    end

    dBdt += -bmort
    d_gain_mort .+= sum(bmort)

    return dBdt, d_gain_mort

end


function zoo_mortality(prms, Z, dZdt, d_gain_mort, t=0)

    zmort = (prms.m_lz .+ prms.m_qz .* Z) .* Z
    dZdt += -zmort
    d_gain_mort .+= sum(zmort)

    return dZdt, d_gain_mort

end


function viral_decay(prms, V, dVdt, d_gain_vde, t=0)

    decay = prms.vde .* V
    dVdt += -decay
    d_gain_vde .+= sum(decay)

    return dVdt, d_gain_vde

end


function total_change_in_d(prms, D, dDdt, d_gain_mort, d_gain_vly, d_gain_vde, t=0)

    dDdt += d_gain_mort .* prms.om_dist_mort
    dDdt += d_gain_vly .* prms.om_dist_lys
    dDdt += d_gain_vde .* prms.om_dist_vde

    dDdt -= D .* prms.rsink

    return dDdt

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

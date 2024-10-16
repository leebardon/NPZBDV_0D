using Printf
using Dates
using NCDatasets
using SparseArrays, LinearAlgebra


function run_CN_test(prms)

    start_time = now()
    trec = prms.nt ÷ prms.nrec

    # Generate empty arrays
    nrec1 = Int(prms.nrec + 1)
    track_n = Array{Float64,2}(undef, prms.nn, nrec1)
    track_c = Array{Float64,2}(undef, prms.nc, nrec1)
    track_b = Array{Float64,2}(undef, prms.nb, nrec1)
    track_don = Array{Float64,2}(undef, prms.ndon, nrec1)
    track_doc = Array{Float64,2}(undef, prms.ndoc, nrec1)
    track_time = Array{Float64,1}(undef, nrec1)

    track_n[:, 1] .= prms.nIC
    track_c[:, 1] .= prms.cIC
    track_b[:, 1] .= prms.bIC
    track_don[:, 1] .= prms.donIC
    track_doc[:, 1] .= prms.docIC

    track_time[1] = 0

    #-------------------- Initial conditions at time t = 0
    btemp = copy(prms.bIC)
    ntemp = copy(prms.nIC)
    ctemp = copy(prms.cIC)
    dontemp = copy(prms.donIC)
    doctemp = copy(prms.docIC)


    for t = 1:prms.nt

        ntemp, ctemp, btemp, dontemp, doctemp = rk4_integrate(ntemp, ctemp, btemp, dontemp, doctemp, prms)

        if mod(t, trec) == 0

            track_n, track_c, track_b, track_don, track_doc, track_time = test_update_tracking_arrs(track_n, track_c, track_b, track_don, track_doc, track_time,
                ntemp, ctemp, btemp, dontemp, doctemp, t, trec, prms)

            tot_n = sum(ntemp) + sum(dontemp) + sum(btemp)
            tot_c = sum(ctemp) + sum(doctemp) + (sum(btemp) * prms.CNr)
            println("Total N: ", tot_n)
            println("Total C: ", tot_c)

        end


        if t == nt
            end_time = now()
            # f = plot(dontemp, doctemp)
            # savefig(f, "dum_fig.jpg")
        end
    end


    return ntemp, ctemp, btemp, dontemp, doctemp, track_time

end


function model_functions(N, C, B, DON, DOC, prms)

    dNdt = zeros(Float64, prms.nn)
    dCdt = zeros(Float64, prms.nc)
    dBdt = zeros(Float64, prms.nb)
    dDONdt = zeros(Float64, prms.ndon)
    dDOCdt = zeros(Float64, prms.ndoc)

    DON_gain_mort = zeros(Float64, 1)
    DOC_gain_mort = zeros(Float64, 1)

    # DIN and DIC supply (turn off rsource and rsink to check if conserving)
    dNdt .+= prms.rsource * (1 / prms.CNr)
    dCdt .+= prms.rsource * (1 - (1 / prms.CNr))

    # bacteria uptake
    dDONdt, dDOCdt, dBdt, dNdt, dCdt = bacteria_uptake(prms, N, B, DON, DOC, dDONdt, dDOCdt, dBdt, dNdt, dCdt)

    #bacterial mortality
    dBdt, DON_gain_mort, DOC_gain_mort = bacterial_mortality(prms, B, dBdt, DON_gain_mort, DOC_gain_mort)

    #sinking rate and split accumulated OM into nd pools
    dDONdt, dDOCdt = total_change_in_d(prms, DON, DOC, dDONdt, dDOCdt, DON_gain_mort, DOC_gain_mort)

    return dNdt, dCdt, dBdt, dDONdt, dDOCdt

end


function bacteria_uptake(prms, N, B, DON, DOC, dDONdt, dDOCdt, dBdt, dNdt, dCdt)
    #NOTE for every 1 N excreted, CNr * C is respired
    # uptake rate controlled by limiting substrate (DOC or DON)

    II, JJ = get_nonzero_axes(prms.CM)

    for j = axes(II, 1)
        # DON and DOC uptake and remineralization (1 - yield)
        muDON = prms.umax_ij[II[j], JJ[j]] .* DON ./ (DON .+ prms.Km_ij[II[j], JJ[j]])
        muDOC = prms.umax_ij[II[j], JJ[j]] .* DOC ./ (DOC .+ prms.Km_ij[II[j], JJ[j]])
        muN = prms.umax_ij[II[j], JJ[j]] .* N ./ (N .+ prms.Km_ij[II[j], JJ[j]])

        # Check if growth limited by DOC or DON+N (what about DIC?)
        yield = prms.y_ij[II[j], JJ[j]]
        mu = min.(yield .* muDOC, muDON .+ muN)
        growth = B[JJ[j], :] .* mu

        mu_Norg = min.(yield .* muDOC, muDON)

        R_NC = 1 / 5  # ratio of N to C in biomass
        uptakeN = min.(muN, max.((mu - mu_Norg), 0)) .* B[JJ[j], :] .* R_NC

        # OM-limited growth rates needed for uptake of OM calc
        growth_Norg = mu_Norg .* B[JJ[j], :] .* R_NC

        EPS = 0.1 #     DUMMY VALUE - to prevent nuerical error for fortran code 
        # calculate the ratios of DON/DOC etc. uptake
        Rup_NC = min.(R_NC, (muDON ./ (muDOC .+ EPS)))

        # Here is uptake and excretion of inorganic nutrients (i.e. sinks in equations for DOC and DON)
        uptakeDOC = growth ./ yield
        uptakeDON = uptakeDOC .* Rup_NC

        # "respDOC" and "respDON" are the remin. sources so go into DIC and NH4 pools 
        respDOC = growth * (1 ./ yield - 1)
        respDON = uptakeDON - growth_Norg

        dBdt[JJ[j], :] .+= growth
        dDONdt[II[j], :] -= uptakeDON
        dDOCdt[II[j], :] -= uptakeDOC

        # proportion remineralized
        dNdt += respDON
        dCdt += respDOC
    end

    return dDONdt, dDOCdt, dBdt, dNdt, dCdt

end

function bacterial_mortality(prms, B, dBdt, DON_gain_mort, DOC_gain_mort)

    bmort = (prms.m_lb .+ prms.m_qb .* B) .* B
    dBdt -= bmort

    DON_gain_mort .+= sum(bmort) * (1 / prms.CNr)
    DOC_gain_mort .+= sum(bmort) * (1 - (1 / prms.CNr))

    return dBdt, DON_gain_mort, DOC_gain_mort

end


function total_change_in_d(prms, DON, DOC, dDONdt, dDOCdt, DON_gain_mort, DOC_gain_mort)

    dDONdt += DON_gain_mort .* prms.om_dist_mort
    dDONdt -= DON .* prms.rsink

    dDOCdt += DOC_gain_mort .* prms.om_dist_mort
    dDOCdt -= DOC .* prms.rsink

    return dDONdt, dDOCdt

end


function get_nonzero_axes(M)

    Cs = sparse(M)
    (II, JJ, _) = findnz(Cs)

    return II, JJ

end




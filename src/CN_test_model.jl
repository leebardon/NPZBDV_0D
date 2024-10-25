using Printf
using Dates, Plots
using NCDatasets
using SparseArrays, LinearAlgebra

global DON_source = 0.06
global DOC_source = 1.0
global OM_sink = 1.0

function run_CN_test(prms)

    start_time = now()
    trec = prms.nt ÷ prms.nrec

    # Generate empty tracking arrays
    nrec1 = Int(prms.nrec + 1)
    trk_b = Array{Float64,2}(undef, prms.nb, nrec1)
    trk_n = Array{Float64,2}(undef, prms.nn, nrec1)
    trk_c = Array{Float64,2}(undef, prms.nc, nrec1)
    trk_don = Array{Float64,2}(undef, prms.ndon, nrec1)
    trk_doc = Array{Float64,2}(undef, prms.ndoc, nrec1)
    trk_c2n = Array{Float64,2}(undef, prms.nn, nrec1)
    trk_time = Array{Float64,1}(undef, nrec1)

    trk_b[:, 1] .= prms.bIC
    trk_n[:, 1] .= prms.nIC
    trk_c[:, 1] .= prms.cIC
    trk_don[:, 1] .= prms.donIC
    trk_doc[:, 1] .= prms.docIC
    trk_c2n[:, 1] .= prms.c2nIC

    trk_time[1] = 0

    #-------------------- Initial conditions at time t = 0
    ntemp = copy(prms.nIC)
    ctemp = copy(prms.cIC)
    btemp = copy(prms.bIC)
    dontemp = copy(prms.donIC)
    doctemp = copy(prms.docIC)

    start_n, start_c = [], []
    # global clim_count = 0
    # global nlim_count = 0

    for t = 1:prms.nt

        ntemp, ctemp, btemp, dontemp, doctemp = rk4_integrate(ntemp, ctemp, btemp, dontemp, doctemp, prms, t)

        if t == 100
            push!(start_n, (sum(ntemp) + sum(dontemp) + sum(btemp)))
            push!(start_c, (sum(ctemp) + sum(doctemp) + sum(btemp) * prms.CNr))
        end

        if mod(t, trec) == 0
            trk_n, trk_c, trk_b, trk_don, trk_doc, trk_c2n, trk_time = test_update_tracking_arrs(trk_n, trk_c, trk_b, trk_don, trk_doc, trk_c2n, trk_time,
                ntemp, ctemp, btemp, dontemp, doctemp, t, trec, prms)

            tot_n = sum(ntemp) + sum(dontemp) + sum(btemp)
            tot_c = sum(ctemp) + sum(doctemp) + (sum(btemp) * prms.CNr)
            println("Total N: ", tot_n)
            println("Total C: ", tot_c)
        end

        if t == nt
            end_time = now()
            println("\nStart N: ", start_n[1])
            println("Start C: ", start_c[1])
            # println("Clim = ", clim_count, "\nNlim = ", nlim_count)
            
            f = plot(trk_time[:], trk_c2n[:], grid=false, label=false, lw=3, xlabel="Days", ylabel="DOC/DON", title="DOC/DON Source = $DOC_source/$DON_source ")
            savefig(f, "dum_fig.png")

            println(trk_c2n[1])
            println(trk_c2n[end])
        end
    end


    return ntemp, ctemp, btemp, dontemp, doctemp, trk_time

end


function model_functions(N, C, B, DON, DOC, prms, t)

    dNdt = zeros(Float64, prms.nn)
    dCdt = zeros(Float64, prms.nc)
    dBdt = zeros(Float64, prms.nb)
    dDONdt = zeros(Float64, prms.ndon)
    dDOCdt = zeros(Float64, prms.ndoc)

    DON_gain_mort = zeros(Float64, 1)
    DOC_gain_mort = zeros(Float64, 1)

    # DIN and DIC supply and sink (set rsource and rsink to zero check if conserving)
    dNdt .+= prms.rsource
    dCdt .+= prms.rsource

    dDONdt .+= DON_source
    dDOCdt .+= DOC_source

    # bacteria uptake
    dDONdt, dDOCdt, dBdt, dNdt, dCdt = bacteria_uptake(prms, N, B, DON, DOC, dDONdt, dDOCdt, dBdt, dNdt, dCdt, t)

    #bacterial mortality
    dBdt, DON_gain_mort, DOC_gain_mort = bacterial_mortality(prms, B, dBdt, DON_gain_mort, DOC_gain_mort)

    #sinking rate and split accumulated OM into nd pools
    dDONdt, dDOCdt = total_change_in_d(prms, DON, DOC, dDONdt, dDOCdt, DON_gain_mort, DOC_gain_mort)

    return dNdt, dCdt, dBdt, dDONdt, dDOCdt

end


function bacteria_uptake(prms, N, B, DON, DOC, dDONdt, dDOCdt, dBdt, dNdt, dCdt, t)

    II, JJ = get_nonzero_axes(prms.CM)

    for j = axes(II, 1)

        yield = prms.y_ij[II[j], JJ[j]]
        umax = prms.umax_ij[II[j], JJ[j]]
        Km = prms.Km_ij[II[j], JJ[j]]
        CNr = prms.CNr

        # Potential uptake rates for DON, DOC and N (DIN) - assume C (DIC) is not limiting
        muDON = umax .* DON ./ (DON .+ Km)
        muDOC = umax .* DOC ./ (DOC .+ Km)

        # Actual growth rate (mu) and production (mu*B) according to limiting resource
        # mu = min.(yield .* muDOC, muDON)

        if (yield .* muDOC) < muDON
            mu = yield .* muDOC
            # global clim_count += 1
        else
            mu = muDON
            # global nlim_count += 1
        end

        growth = B[JJ[j], :] .* mu

        # calculate whether the uptake ratio of N:C is limited by NC of biomass (1/CNr) of env (muDON/muDOC)
        up_NC = min.(1 / CNr, (muDON ./ muDOC)) # should probs use up_CN for clarity)

        # Uptake and excretion (i.e. sinks in equations for DOC and DON)
        uptakeDOC = (growth .* (1 ./ up_NC)) ./ yield
        uptakeDON = uptakeDOC .* up_NC
        respDOC = (growth .* (1 ./ up_NC)) .* (1 ./ yield - 1)
        respDON = uptakeDON - growth

        dBdt[JJ[j], :] += growth
        dDONdt[II[j], :] -= uptakeDON
        dDOCdt[II[j], :] -= uptakeDOC
        dNdt += respDON
        dCdt += respDOC
    end

    return dDONdt, dDOCdt, dBdt, dNdt, dCdt

end

function bacterial_mortality(prms, B, dBdt, DON_gain_mort, DOC_gain_mort)

    bmort = (prms.m_lb .+ prms.m_qb .* B) .* B
    dBdt -= bmort

    DON_gain_mort .+= sum(bmort)
    DOC_gain_mort .+= sum(bmort) .* prms.CNr

    return dBdt, DON_gain_mort, DOC_gain_mort

end


function total_change_in_d(prms, DON, DOC, dDONdt, dDOCdt, DON_gain_mort, DOC_gain_mort)

    dDONdt += DON_gain_mort .* prms.om_dist_mort
    # dDONdt -= DON .* prms.rsink
    dDONdt -= DON .* OM_sink

    dDOCdt += DOC_gain_mort .* prms.om_dist_mort
    # dDOCdt -= DOC .* prms.rsink
    dDOCdt -= DOC .* OM_sink

    return dDONdt, dDOCdt

end


function get_nonzero_axes(M)

    Cs = sparse(M)
    (II, JJ, _) = findnz(Cs)

    return II, JJ

end



# For DIN
# muN = umax .* N ./ (N .+ Km)
# mu = min.(yield .* muDOC, muDON .+ muN)

# OM-limited growth rates needed for uptake of OM calc
# mu_Norg = min.(yield .* muDOC, muDON)
# growth_Norg = B[JJ[j], :] .* mu_Norg
# uptakeN = min.(muN, max.((mu - mu_OM), 0)) .* B[JJ[j], :] .* (1/CNr)
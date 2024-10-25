
include("CN_test_model.jl")

function test_update_tracking_arrs(track_n, track_c, track_b, track_don, track_doc, track_c2n, track_time, ntemp, ctemp, btemp, dontemp, doctemp, t, trec, prms)

    j = Int(t ÷ trec + 1)
    t_id = t .* prms.dt
    track_b[:, j] .= btemp
    track_n[:, j] .= ntemp
    track_c[:, j] .= ctemp
    track_don[:, j] .= dontemp
    track_doc[:, j] .= doctemp
    track_c2n[:, j] .= doctemp ./ dontemp
    track_time[j] = t_id

    @printf("Day %7.1f out of %5.0f = %4.0f%% done at %s \n", t_id, prms.days, t_id / prms.days * 100, now())

    return track_n, track_c, track_b, track_don, track_doc, track_c2n, track_time

end


function rk4_integrate(ntemp, ctemp, btemp, dontemp, doctemp, prms, t)

    dXdt1 = model_functions(ntemp, ctemp, btemp, dontemp, doctemp, prms, t)

    track_n1 = ntemp .+ prms.dt / 2 .* dXdt1[1]
    track_c1 = ctemp .+ prms.dt / 2 .* dXdt1[2]
    track_b1 = btemp .+ prms.dt / 2 .* dXdt1[3]
    track_dn1 = dontemp .+ prms.dt / 2 .* dXdt1[4]
    track_dc1 = doctemp .+ prms.dt / 2 .* dXdt1[5]

    dXdt2 = model_functions(track_n1, track_c1, track_b1, track_dn1, track_dc1, prms, t)

    track_n2 = ntemp .+ prms.dt / 2 .* dXdt2[1]
    track_c2 = ctemp .+ prms.dt / 2 .* dXdt2[2]
    track_b2 = btemp .+ prms.dt / 2 .* dXdt2[3]
    track_dn2 = dontemp .+ prms.dt / 2 .* dXdt2[4]
    track_dc2 = doctemp .+ prms.dt / 2 .* dXdt2[5]

    dXdt3 = model_functions(track_n2, track_c2, track_b2, track_dn2, track_dc2, prms, t)

    track_n3 = ntemp .+ prms.dt .* dXdt3[1]
    track_c3 = ctemp .+ prms.dt .* dXdt3[2]
    track_b3 = btemp .+ prms.dt .* dXdt3[3]
    track_dn3 = dontemp .+ prms.dt .* dXdt3[4]
    track_dc3 = doctemp .+ prms.dt .* dXdt3[5]

    dXdt4 = model_functions(track_n3, track_c3, track_b3, track_dn3, track_dc3, prms, t)

    ntemp .+= (dXdt1[1] .+ 2 .* dXdt2[1] .+ 2 .* dXdt3[1] .+ dXdt4[1]) .* (prms.dt / 6)
    ctemp .+= (dXdt1[2] .+ 2 .* dXdt2[2] .+ 2 .* dXdt3[2] .+ dXdt4[2]) .* (prms.dt / 6)
    btemp .+= (dXdt1[3] .+ 2 .* dXdt2[3] .+ 2 .* dXdt3[3] .+ dXdt4[3]) .* (prms.dt / 6)
    dontemp .+= (dXdt1[4] .+ 2 .* dXdt2[4] .+ 2 .* dXdt3[4] .+ dXdt4[4]) .* (prms.dt / 6)
    doctemp .+= (dXdt1[5] .+ 2 .* dXdt2[5] .+ 2 .* dXdt3[5] .+ dXdt4[5]) .* (prms.dt / 6)


    return ntemp, ctemp, btemp, dontemp, doctemp


end

###################################################
#SUBSTRATE TRAITS ("averages")

function ordered_uptake_arr(n, phy=false)

    max_uptake_i = collect(10 .^ range(-2, stop=1, length=n))

    if phy == true
        max_uptake_i = collect(10 .^ range(0, stop=1, length=n))
    end

    return max_uptake_i

end


function random_uptake_arr(n)
    # randomly selecting along a log range 
    xhi = log10(1e1)
    xlo = log10(1e-1)
    max_uptake_i = 10 .^ ((rand(n) .- 0.5) .* (xhi - xlo) .+ mean([xhi xlo]))

    return max_uptake_i

end


function reproducible_Fg(n; rng=GlobalRNG)

    return rand(rng, n)

end


function apply_tradeoff(nconsumers, nresources, CM, max_i, season, run_type)
    #fraction of the proteome devoted to growth (F_g) vs. affinity (F_a)
    #NOTE could also explore growth-defense tradeoff https://www.nature.com/articles/s41396-020-0619-1 or growth-yield
    if run_type == 1
        Fg_b, Fg_p = reproducible_Fg(nconsumers), reproducible_Fg(nconsumers)
        Fa_b, Fa_p = 1.0 .- Fg_b, 1.0 .- Fg_p
    elseif run_type == 2
        Fg_b, Fg_p = get_prescribed_params("Fg_b"), get_prescribed_params("Fg_p")
        Fa_b, Fa_p = 1.0 .- Fg_b, 1.0 .- Fg_p
    end

    if nresources == 1
        vmax_ij = set_vmax_ij(nresources, nconsumers, max_i, Fg_p)
        Kp_ij = set_Kp_ij(nresources, nconsumers, Fa_p, CM, vmax_ij)
        return vmax_ij, Kp_ij, Fg_p
    else
        umax_ij = set_umax_ij(nresources, nconsumers, max_i, CM, Fg_b)
        Km_ij = set_Km_ij(nresources, nconsumers, Fa_b, CM, umax_ij)
        return umax_ij, Km_ij, Fg_b
    end

end


function set_umax_ij(nd, nb, umax_i, CM, F_g)
    # Max growth of b on d (note - if F_g == 0.5 and F_a is 0.5, then max_ij == max_i)
    umax_ij = zeros(nd, nb)
    for i = 1:nd
        umax_ij[i, :] = (umax_i[i] * F_g) ./ (0.5 .* CM[i, :])
    end

    return umax_ij
end


function set_vmax_ij(nn, np, vmax_i, F_g)
    # Max growth of p on n (note - only applicable to 1N runs, rewrite like B on D if more N needed)
    vmax_ij = zeros(nn, np)
    for i = 1:nn
        vmax_ij[i, :] = (vmax_i .* F_g) ./ 0.5
    end

    return vmax_ij
end


function set_Km_ij(nd, nb, F_a, CM, umax_ij)
    # half-sat of b on d (note, if F_a is 0.5, then Km_ij == Km_i)
    affinity = zeros(nd, nb)
    Km_ij = zeros(nd, nb)
    for i = 1:nd
        affinity[i, :] = 10 * F_a ./ 0.5 .* CM[i, :]
        Km_ij[i, :] = umax_ij[i, :] ./ affinity[i, :]
    end

    return Km_ij
end


function set_Kp_ij(nn, np, F_a, CM, vmax_ij)
    # half-sat of p on n
    affinity = zeros(nn, np)
    Kp_ij = zeros(nn, np)
    for i = 1:nn
        affinity[i, :] = 10 * F_a ./ 0.5 .* CM[i, :]
        Kp_ij[i, :] = vmax_ij[i, :] ./ affinity[i, :]
    end

    return Kp_ij
end




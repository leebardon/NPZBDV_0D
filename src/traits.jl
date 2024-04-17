###################################################
#SUBSTRATE TRAITS ("averages")

function ordered_uptake_arr(n, phy=false)

    max_uptake_i = collect(10 .^ range(-2,stop=1,length=n))

    if phy==true
        max_uptake_i = collect(10 .^ range(0,stop=1,length=n))
    end

    return max_uptake_i

end


function random_uptake_arr(n)
    # randomly selecting along a log range 
    xhi = log10(1e1)
    xlo = log10(1e-1) 
    max_uptake_i = 10 .^ ((rand(n).-0.5).*(xhi-xlo).+mean([xhi xlo]))

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
        Fa_b, Fa_p = 1. .- Fg_b, 1. .- Fg_p
    elseif run_type == 2
        Fg_b, Fg_p = get_prescribed_params("Fg_b"), get_prescribed_params("Fg_p") 
        Fa_b, Fa_p = 1. .- Fg_b, 1. .- Fg_p
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
        umax_ij[i,:] = (umax_i[i] * F_g) ./ (0.5 .* CM[i,:])  
    end

    return umax_ij
end 


function set_vmax_ij(nn, np, vmax_i, F_g)
    # Max growth of p on n (note - only applicable to 1N runs, rewrite like B on D if more N needed)
    vmax_ij = zeros(nn, np)
    for i = 1:nn
        vmax_ij[i,:] = (vmax_i .* F_g) ./ 0.5 
    end

    return vmax_ij
end 


function set_Km_ij(nd, nb, F_a, CM, umax_ij)
    # half-sat of b on d (note, if F_a is 0.5, then Km_ij == Km_i)
    affinity = zeros(nd, nb)
    Km_ij = zeros(nd, nb)
        for i = 1:nd
        affinity[i,:] = 10*F_a ./ 0.5.*CM[i,:]
        Km_ij[i,:] = umax_ij[i,:] ./ affinity[i,:] 
    end

    return Km_ij
end


function set_Kp_ij(nn, np, F_a, CM, vmax_ij)
    # half-sat of p on n
    affinity = zeros(nn, np)
    Kp_ij = zeros(nn, np)
    for i = 1:nn
        affinity[i,:] = 10*F_a ./ 0.5.*CM[i,:]
        Kp_ij[i,:] = vmax_ij[i,:] ./ affinity[i,:] 
    end

    return Kp_ij
end




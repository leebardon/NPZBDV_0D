

function rk4(ntemp, ptemp, ztemp, btemp, dtemp, vtemp, prms, t, lysis)

    n, p, z, b, d, v = copy(ntemp), copy(ptemp), copy(ztemp), copy(btemp), copy(dtemp), copy(vtemp)

    dXdt1 = model_functions(n, p, z, b, d, v, prms, t, lysis)
    
    track_n1 = n .+ prms.dt/2 .* dXdt1[1]
    track_p1 = p .+ prms.dt/2 .* dXdt1[2]
    track_z1 = z .+ prms.dt/2 .* dXdt1[3]
    track_b1 = b .+ prms.dt/2 .* dXdt1[4]
    track_d1 = d .+ prms.dt/2 .* dXdt1[5]
    track_v1 = v .+ prms.dt/2 .* dXdt1[6]

    dNdt2, dPdt2, dZdt2, dBdt2, dDdt2, dVdt2 = model_functions(track_n1, track_p1, track_z1, track_b1, track_d1, track_v1, prms, t, lysis)

    track_n2 = n .+ prms.dt/2 .* dNdt2
    track_p2 = p .+ prms.dt/2 .* dPdt2
    track_z2 = z .+ prms.dt/2 .* dZdt2
    track_b2 = b .+ prms.dt/2 .* dBdt2
    track_d2 = d .+ prms.dt/2 .* dDdt2
    track_v2 = v .+ prms.dt/2 .* dVdt2
    
    dNdt3, dPdt3, dZdt3, dBdt3, dDdt3, dVdt3  = model_functions(track_n2, track_p2, track_z2, track_b2, track_d2, track_v2, prms, t, lysis)

    track_n3 = n .+ prms.dt .* dNdt3
    track_p3 = p .+ prms.dt .* dPdt3
    track_z3 = z .+ prms.dt .* dZdt3
    track_b3 = b .+ prms.dt .* dBdt3
    track_d3 = d .+ prms.dt .* dDdt3
    track_v3 = v .+ prms.dt .* dVdt3

    dNdt4, dPdt4, dZdt4, dBdt4, dDdt4, dVdt4 = model_functions(track_n3, track_p3, track_z3, track_b3, track_d3, track_v3, prms, t, lysis)

    n .+= (dXdt1[1] .+ 2 .* dNdt2 .+ 2 .* dNdt3 .+ dNdt4) .* (prms.dt / 6)
    p .+= (dXdt1[2] .+ 2 .* dPdt2 .+ 2 .* dPdt3 .+ dPdt4) .* (prms.dt / 6)
    z .+= (dXdt1[3] .+ 2 .* dZdt2 .+ 2 .* dZdt3 .+ dZdt4) .* (prms.dt / 6)
    b .+= (dXdt1[4].+ 2 .* dBdt2 .+ 2 .* dBdt3 .+ dBdt4) .* (prms.dt / 6)
    d .+= (dXdt1[5] .+ 2 .* dDdt2 .+ 2 .* dDdt3 .+ dDdt4) .* (prms.dt / 6)
    v .+= (dXdt1[6] .+ 2 .* dVdt2 .+ 2 .* dVdt3 .+ dVdt4) .* (prms.dt / 6)


    return n, p, z, b, d, v

end
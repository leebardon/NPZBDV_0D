

function rk4(ntemp, ctemp, ptemp, ztemp, btemp, dtemp, vtemp, prms, t, lysis)


    dXdt1 = model_functions(ntemp, ctemp, ptemp, ztemp, btemp, dtemp, vtemp, prms, t, lysis)

    track_n1 = ntemp .+ prms.dt/2 .* dXdt1[1]
    track_c1 = ctemp .+ prms.dt/2 .* dXdt1[2]
    track_p1 = ptemp .+ prms.dt/2 .* dXdt1[3]
    track_z1 = ztemp .+ prms.dt/2 .* dXdt1[4]
    track_b1 = btemp .+ prms.dt/2 .* dXdt1[5]
    track_d1 = dtemp .+ prms.dt/2 .* dXdt1[6]
    track_v1 = vtemp .+ prms.dt/2 .* dXdt1[7]

    # N1tot = sum(track_p1) + sum(track_b1) + sum(track_z1) + sum(track_n1) + sum(track_d1)
    # nan_or_inf(N1tot) && @error("Nan or Inf at t=$t: \n N1tot: $N1tot")
    
    dXdt2 = model_functions(track_n1, track_c1, track_p1, track_z1, track_b1, track_d1, track_v1, prms, t, lysis)

    track_n2 = ntemp .+ prms.dt/2 .* dXdt2[1]
    track_c2 = ctemp .+ prms.dt/2 .* dXdt2[2]
    track_p2 = ptemp .+ prms.dt/2 .* dXdt2[3]
    track_z2 = ztemp .+ prms.dt/2 .* dXdt2[4]
    track_b2 = btemp .+ prms.dt/2 .* dXdt2[5]
    track_d2 = dtemp .+ prms.dt/2 .* dXdt2[6]
    track_v2 = vtemp .+ prms.dt/2 .* dXdt2[7]

    # N2tot = sum(track_p2) + sum(track_b2) + sum(track_z2) + sum(track_n2) + sum(track_d2)
    # nan_or_inf(N2tot) && @error("Nan or Inf at t=$t: \n N2tot: $N2tot") 
    
    dXdt3 = model_functions(track_n2, track_c2, track_p2, track_z2, track_b2, track_d2, track_v2, prms, t, lysis)

    track_n3 = ntemp .+ prms.dt .* dXdt3[1]
    track_c3 = ctemp .+ prms.dt .* dXdt3[2]
    track_p3 = ptemp .+ prms.dt .* dXdt3[3]
    track_z3 = ztemp .+ prms.dt .* dXdt3[4]
    track_b3 = btemp .+ prms.dt .* dXdt3[5]
    track_d3 = dtemp .+ prms.dt .* dXdt3[6]
    track_v3 = vtemp .+ prms.dt .* dXdt3[7]

    # N3tot = sum(track_p3) + sum(track_b3) + sum(track_z3) + sum(track_n3) + sum(track_d3)
    # nan_or_inf(N3tot) && @error("Nan or Inf at t=$t: \n N3tot: $N3tot") 

    dXdt4 = model_functions(track_n3, track_c3, track_p3, track_z3, track_b3, track_d3, track_v3, prms, t, lysis)

    ntemp .+= (dXdt1[1] .+ 2 .* dXdt2[1] .+ 2 .* dXdt3[1] .+ dXdt4[1]) .* (prms.dt / 6)
    ctemp .+= (dXdt1[2] .+ 2 .* dXdt2[2] .+ 2 .* dXdt3[2] .+ dXdt4[2]) .* (prms.dt / 6)
    ptemp .+= (dXdt1[3] .+ 2 .* dXdt2[3] .+ 2 .* dXdt3[3] .+ dXdt4[3]) .* (prms.dt / 6)
    ztemp .+= (dXdt1[4] .+ 2 .* dXdt2[4] .+ 2 .* dXdt3[4] .+ dXdt4[4]) .* (prms.dt / 6)
    btemp .+= (dXdt1[5] .+ 2 .* dXdt2[5] .+ 2 .* dXdt3[5] .+ dXdt4[5]) .* (prms.dt / 6)
    dtemp .+= (dXdt1[6] .+ 2 .* dXdt2[6] .+ 2 .* dXdt3[6] .+ dXdt4[6]) .* (prms.dt / 6)
    vtemp .+= (dXdt1[7] .+ 2 .* dXdt2[7] .+ 2 .* dXdt3[7] .+ dXdt4[7]) .* (prms.dt / 6)
     
    # Ntot = sum(ptemp) + sum(btemp) + sum(ntemp) + sum(dtemp) + sum(ztemp) + sum(vtemp)
    # nan_or_inf(Ntot) && @error("Nan or Inf at t=$t: \n Ntot: $Ntot") 

    return ntemp, ctemp, ptemp, ztemp, btemp, dtemp, vtemp


end
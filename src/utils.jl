
import REPL
using REPL.TerminalMenus
using Printf, NCDatasets
using Dates
using SparseArrays



function message(v::String, nd::Int64=0, nb::Int64=0, np::Int64=0, nz::Int64=0, fsaven::String="")

    m = Dict(
        "LY1" => "\nSelect viral lysis type:",
        "LY2" => ["Explicit", "Implicit"],
        "GZ1" => "\nInclude grazers?",
        "GZ2" => ["yes", "no"],
        "CB1" => "\nInclude carbon pools?",
        "CB2" => ["yes", "no"],
        "TM" => "\nEnter simulation runtime (years):",
    )

    return m["$v"]

end


function set_savefiles(launch_time, years, nn, nc, np, nz, nb, ndn, ndc, nv, lysis=0)

    fsave = "results/outfiles/"

    if lysis != 1
        fsaven = string(fsave, Dates.format(launch_time, "yymmdd_HH:MM"), "_$(years)y_$(nn)N$(nc)C$(np)P$(nz)Z$(nb)B$(ndn+ndc)D$(nv)V.nc")
    else
        fsaven = string(fsave, Dates.format(launch_time, "yymmdd_HH:MM"), "_$(years)y_$(nn)N$(nc)C$(np)P$(nz)Z$(nb)B$(ndn+ndc)D$(nv)V_L.nc")
    end

    return fsaven

end


function set_logger(launch_time)

    loginfo = string(Dates.format(launch_time, "yyyymmdd_HHMM"), ".log")
    logger = activate_logger(loginfo)

    return logger

end


function activate_logger(loginfo)

    logger = TeeLogger(
        MinLevelLogger(FileLogger("logs/$loginfo"), Logging.Info),
        MinLevelLogger(FileLogger("logs/error.log"), Logging.Warn),
    ); 
    
    global_logger(logger)

    return logger 

end

function log_params(prms, lysis)

    @info(
    """Model Params: 
    days:           $(prms.days)
    ts/day (dt):    $(prms.dt)
    total ts (nt):  $(prms.nt)                    
    nn:             $(prms.nn)                    
    np:             $(prms.np)                     
    nz:             $(prms.nz)        
    nb:             $(prms.nb)         
    ndn:            $(prms.ndn) 
    ndc:            $(prms.ndc) 
    nv:             $(prms.nv)   
    nIC:            $(prms.nIC[1])                        
    pIC:            $(prms.pIC[1])
    zIC:            $(prms.zIC[1])
    bIC:            $(prms.bIC[1])
    dnIC:           $(prms.dnIC[1])
    dcIC:           $(prms.dcIC[1])
    vIC:            $(prms.vIC[1])
    vmax_i:         $(prms.vmax_i)
    vmax_ij:        $(prms.vmax_ij)
    umax_i:         $(prms.umax_i)
    umax_ij:        $(prms.umax_ij)
    Kp_i:           $(prms.Kp_i)
    Kp_ij:          $(prms.Kp_ij)
    Km_i:           $(prms.Km_i)
    Km_ij:          $(prms.Km_ij)
    y_ij:           $(prms.y_ij[1])
    K_g:            $(prms.K_g[1])
    g_max:          $(prms.g_max[1])
    γ:              $(prms.γ[1])
    m_lp:           $(prms.m_lp[1])
    m_qp:           $(prms.m_qp[1])
    m_lb:           $(prms.m_lb[1])
    m_qb:           $(prms.m_qb[1])
    m_lz:           $(prms.m_lz[1])
    m_qz:           $(prms.m_qz[1])
    rsource:        $(prms.rsource)
    rsink:          $(prms.rsink)
    om_dist_mort:   $(prms.om_dist_mort)  
    om_dist_lys:    $(prms.om_dist_lys)   
    om_dist_vde:    $(prms.om_dist_vde) 
    Fg_p:           $(prms.Fg_p)
    Fg_b:           $(prms.Fg_b) 
    lysis:          $(lysis)     # 1 = explicit, 2 = implicit
    savefile:       $(prms.fsaven)
    """) 
    
    print_info(prms)

end


function print_info(prms)

    @printf("\n np = %5.0f \n nb = %5.0f \n nz = %5.0f \n nn = %5.0f \n ndn = %5.0f \n ndc = %5.0f \n nv = %5.0f \n days = %5.0f \n\n", prms.np, prms.nb, prms.nz, prms.nn, prms.ndn, prms.ndc, prms.nv, prms.days)
    println("File will be saved as: ", prms.fsaven)
    println("nt = ", prms.nt)

end



function bacteria_num()
    println(message("BA"))
    input = readline()
    nb = parse(Int64, input)
    return nb
end


function user_select()

    println(message("T"))
    input = readline()
    tt = parse(Int64, input) 

    println(message("REC"))
    input = readline()
    nrec = parse(Int64, input) 

    println(message("OM"))
    input = readline()
    nd = parse(Int64, input) 

    nb = bacteria_num()
    if !iseven(nb)
        println("\n !!! PLEASE SELECT EVEN NUMBER OF BACTERIA !!! \n\n")
        nb = bacteria_num()
    end

    println(message("NPZ"))
    np = 1
    nz = 2
    nn = 1

    yield = request(message("Y2"), RadioMenu(message("Y1")))
    OM_supply_weight = request(message("SW2"), RadioMenu(message("SW1")))
    println(message("SUB"))
    uptake = request(message("UP2"), RadioMenu(message("UP1")))

    #TODO add error handling to catch int in input

    println(message("ENV"))
    println(message("S1"))
    input = readline()
    rsource = parse(Float64, input) 

    println(message("S2"))
    input = readline()
    rsink = parse(Float64, input) 
    
    pulse = request(message("PU2"), RadioMenu(message("PU1")))

    return tt, nrec, nd, nb, np, nz, nn, yield, OM_supply_weight, uptake, rsource, rsink, pulse

end


function load_matrix(mtype, nd, nb, np=0, nz=0)
    
    M = jldopen("results/matrices/$(mtype)_$(np)p$(nd)d$(nb)b$(nz)z.jdl", "r") do file
        read(file, "A")
    end

    return M

end


function save_matrices(M1, M2, nd, nb, np, nz)

    jldopen("results/matrices/CM_$(np)p$(nd)d$(nb)b$(nz)z.jdl", "w") do file
        write(file, "A", M1)  
    end
    jldopen("results/matrices/GrM_$(np)p$(nd)d$(nb)b$(nz)z.jdl", "w") do file
            write(file, "A", M2)  
    end

end


function test_vals(arr)

    e_msg = "\n Nan or inf found in timestep "
    w_msg = "\n Weird values found in timestep "
    check = run_checks(arr)

    return check, arr, e_msg, w_msg

end


function run_checks(vals)

    for x in vals
        if nan_or_inf(x)
            return "e"
        elseif big_or_small(x)
            return "w"
        end
    end

end


function nan_or_inf(x)

    if typeof(x) == Float64 || typeof(x) == Int64
        if isnan(x) || !isfinite(x)
            return true
        end
    else
        if any(isnan.(x)) || any(isinf.(x))
            return true
        end
    end

    return false

end

function big_or_small(x)

    if typeof(x) == Float64 || typeof(x) == Int64
        if x > 1e10 || x < -1e10 
            return true
        end
    else
        for i in x
            if i > 1e10 || i < -1e10
                return true
            end
        end
    end

    return false

end


function get_endpoints(vars, ds=nothing)

    endpoints = Vector{Any}()

    if ds !== nothing
        for v in vars
            append!(endpoints, [ds["$v"][:,end]])
        end
    else
        for v in vars
            append!(endpoints, [v[:,end]])
        end
    end

    return endpoints
end



function print_info(start_time, prms, nt)

    @printf("\n np = %5.0f \n nb = %5.0f \n nz = %5.0f \n nn = %5.0f \n nd = %5.0f \n T = %5.0f \n\n", prms.np, prms.nb, prms.nz, prms.nn,  prms.nd, prms.tt)
    open("jlog.txt","w") do f
        write(f,@sprintf("np = %5.0f, nb = %5.0f, nz = %5.0f, nn = %5.0f, nd = %5.0f, T = %5.0f \n", prms.np, prms.nb, prms.nz, prms.nn,  prms.nd, prms.tt))
    end

    fsaven = string(prms.fsave,"_", Dates.format(start_time, "yyyymmdd"), ".nc")
    if isfile(fsaven)
        fsaven = string(prms.fsave, "_", Dates.format(start_time, "yyyymmdd_HHMM"), ".nc")
    end

    println("Starting time: $start_time, \nFile will be saved as: $fsaven")

    println("\nConsumption Matrix (CM):")
    display("text/plain", prms.CM)

    println("\nGrazing Matrix (GrM):")
    display("text/plain", prms.GrM)

    println("nt = ", nt)

    return fsaven

end


function update_tracking_arrs(track_n, track_c, track_p, track_z, track_b, track_dn, track_dc, track_v, track_time, ntemp, ctemp, ptemp, ztemp, btemp, dntemp, dctemp, vtemp, t, trec, prms)

    j = Int(t÷trec + 1)
    t_id = t.*prms.dt
    track_p[:,j] .= ptemp
    track_b[:,j] .= btemp 
    track_z[:,j] .= ztemp 
    track_n[:,j] .= ntemp 
    track_c[:,j] .= ctemp 
    track_dn[:,j] .= dntemp
    track_dc[:,j] .= dctemp
    track_v[:,j] .= vtemp
    track_time[j] = t_id 

    @printf("Day %7.1f out of %5.0f = %4.0f%% done at %s \n", t_id, prms.days, t_id/prms.days*100, now())

    return track_n, track_c, track_p, track_z, track_b, track_dn, track_dc, track_v, track_time

end


function savetoNC(fsaven, p, b, z, n, c, dn, dc, v, timet, tst, tfn, prms, pulse)

    outdir = "/home/lee/Dropbox/Development/NPZBDV_0D/"
    path = joinpath(outdir, fsaven) 
    println("Saving to: ", path)

    f = NCDataset(path, "c") #c for create

    # define the dim of p, b, z, n, dn, dc
    defDim(f,"np", prms.np)
    defDim(f,"nb", prms.nb)
    defDim(f,"nz", prms.nz)
    defDim(f,"nn", prms.nn)
    defDim(f,"nc", prms.nc)
    defDim(f,"ndn", prms.ndn)
    defDim(f,"ndc", prms.ndc)
    defDim(f,"nv", prms.nv)

    # define the dim of the time length
    nrec1 = Int(prms.nrec+1) #bc i added time 0
    nprey = prms.np + prms.nb
    
    defDim(f,"nrec",nrec1)
    defDim(f,"nprey",nprey)
   
    # info
    f.attrib["title"] = "NPZBDV 0D model i/o"
    f.attrib["Start time"] = string(tst)
    f.attrib["End time"] = string(tfn)
    f.attrib["Run time"] = string(tfn - tst) 

     # simulated results
     w = defVar(f,"p",Float64,("np","nrec"))
     w[:,:] = p
     w.attrib["units"] = "mmol/m3 N biomass"
 
     w = defVar(f,"b",Float64,("nb","nrec"))
     w[:,:] = b
     w.attrib["units"] = "mmol/m3 N biomass"
 
     w = defVar(f,"z",Float64,("nz","nrec"))
     w[:,:] = z
     w.attrib["units"] = "mmol/m3 N biomass"
     
     w = defVar(f,"n",Float64,("nn","nrec"))
     w[:,:] = n
     w.attrib["units"] = "mmol/m3 N OM"

     w = defVar(f,"c",Float64,("nc","nrec"))
     w[:,:] = c
     w.attrib["units"] = "mmol/m3 C OM"
 
     w = defVar(f,"dn",Float64,("ndn","nrec"))
     w[:,:] = dn
     w.attrib["units"] = "mmol/m3 N OM"

     w = defVar(f,"dc",Float64,("ndc","nrec"))
     w[:,:] = dc
     w.attrib["units"] = "mmol/m3 C OM"

     w = defVar(f,"v",Float64,("nv","nrec"))
     w[:,:] = v
     w.attrib["units"] = "mmol/m3 N OM"
     
    #  w = defVar(f,"uptake",Float64,("nd","nb"))
    #  w[:,:] = uptake
    #  w.attrib["units"] = "mmol/m3 N per d; uptake matrix"
     
     # w = defVar(f,"pIC",Float64,("np",))
     # w[:,:] = prms.pIC
     # w.attrib["units"] = "mmol/m3 N biomass"
 
     # w = defVar(f,"bIC",Float64,("nb",))
     # w[:,:] = prms.bIC
     # w.attrib["units"] = "mmol/m3 N biomass"
 
     # w = defVar(f,"zIC",Float64,("nz",))
     # w[:,:] = prms.zIC
     # w.attrib["units"] = "mmol/m3 N biomass"
 
     # w = defVar(f,"nIC",Float64,("nn",))
     # w[:,:] = prms.nIC
     # w.attrib["units"] = "mmol/m3 N OM"
 
     # w = defVar(f,"dIC",Float64,("nd",))
     # w[:,:] = prms.dIC
     # w.attrib["units"] = "mmol/m3 N OM"
     
     w = defVar(f, "timet", Float64, ("nrec",))
     w[:] = timet
     w.attrib["units"] = "days"
 
    #  w = defVar(f, "K_n", Float64, ("np",))
    #  w[:] = prms.K_n
    #  w.attrib["units"] = "m3/mmol; half-sat rate of p"
 
     w = defVar(f,"CM",Float64,("ndn","nb"))
     w[:,:] = prms.CM
     w.attrib["units"] = "Consumption Matrix"
 
     w = defVar(f,"GrM",Float64,("nz","nprey"))
     w[:,:] = prms.GrM
     w.attrib["units"] = "Grazing Matrix"
 
     # w = defVar(f, "y_i", Float64, ("nd", "nb"))
     # w[:] = prms.y_i
     # w.attrib["units"] = "mol B/mol C; yield"
 
     # w = defVar(f, "y_ij", Float64, ("nd","nb"))
     # w[:,:] = prms.y_ij
     # w.attrib["units"] = "per d; max yield rate"
 
     # w = defVar(f, "umax_i", Float64, ("nd",))
     # w[:,:] = prms.umax_i
     # w.attrib["units"] = "per d; max uptake rate"
     
     # w = defVar(f, "umax_ij", Float64, ("nd", "nb"))
     # w[:,:] = prms.umax_ij
     # w.attrib["units"] = "per d; max uptake rate"
     
     # w = defVar(f, "Km_ij", Float64, ("nd", "nb"))
     # w[:] = prms.Km_ij
     # w.attrib["units"] = "mmol/m3; half-sat"
     
     # w = defVar(f, "m_lb", Float64, ("nb",))
     # w[:,:] = prms.m_lb
     # w.attrib["units"] = "m3/mmol; linear death rate of b"
 
     # w = defVar(f, "m_qb", Float64, ("nb",))
     # w[:,:] = prms.m_qb
     # w.attrib["units"] = "m3/mmol; quadratic death rate of b"
     
     # w = defVar(f, "K_g", Float64, ("nz",))
     # w[:] = prms.K_g
     # w.attrib["units"] = "m3/mmol; half-sat rate of z"
     
     # w = defVar(f, "γ", Float64, ("nz",))
     # w[:] = prms.γ
     # w.attrib["units"] = "fraction of digestion of z"
     
     # w = defVar(f, "m_lz", Float64, ("nz",))
     # w[:,:] = prms.m_lz
     # w.attrib["units"] = "m3/mmol; linear death rate of z"
 
     # w = defVar(f, "m_qz", Float64, ("nz",))
     # w[:,:] = prms.m_qz
     # w.attrib["units"] = "m3/mmol; quadratic death rate of z"
     
     close(f)
 
 end

# function define_dims(ds, prms, nrec1)

#     vars = Dict("nb" => prms.nb, "nd" => prms.nd, "nrec" => nrec1)

#     for (k, v) in x
#         defDim(ds, k, v)
#     end

#     return ds

# end

# BELOW WAS INSERTED INTO FUNCTIONS TO TRACE NAN AND INF WEIRDNESS - CAUSED BY USING UNDEF TO 
# INITIALISE EMPTY ARRS TO BE LATER USED DURING INTEGRATION


# check, data, e_msg, w_msg = test_vals([dNdt, dPdt, dZdt, dBdt, dDdt])
# if check=="e"
#     @error("$e_msg $t at j=$j: \n $data \n")
# elseif check=="w"
#     print("warn")
#     @error("$w_msg $t at j=$j: \n $data \n")
# end

#         #TODO find a better way to do this so can traceback source of fault and not need to repeat code blocks
#         #NOTE this is probably something that can be improved with proper unit testing
#         check, data, e_msg, w_msg = test_vals([uptake, mort, dNdt, dPdt, d_gain_total])
#         if check=="e"
#             @error("$e_msg $t at i=$i: \n $data \n")
#         elseif check=="w"
#             print("warn")
#             @error("$w_msg $t at i=$i: \n $data \n ")
#         end


#         check, data, e_msg, w_msg = test_vals([uptake, dDdt, dBdt, dNdt])
#         if check=="e"
#             @error("$e_msg $t at j=$j: \n $data \n")
#         elseif check=="w"
#             print("warn")
#             @error("$w_msg $t at j=$j: \n $data \n")
#         end


#         check, data, e_msg, w_msg = test_vals([prey, g, dZdt, dNdt, dPdt])
#         if check=="e"
#             @error("$e_msg $t at k=$k: \n $data \n")
#         elseif check=="w"
#             print("warn")
#             @error("$w_msg $t at k=$k: \n $data \n")
#         end  


#         check, data, e_msg, w_msg = test_vals([prey, g, dZdt, dNdt, dBdt])
#         if check=="e"
#             @error("$e_msg $t at k=$k: \n $data \n")
#         elseif check=="w"
#             print("warn")
#             @error("$w_msg $t at k=$k: \n $data \n")
#         end


#         check, data, e_msg, w_msg = test_vals([bmort, dBdt, d_gain_total])
#         if check=="e"
#             @error("$e_msg $t : \n $data \n")
#         elseif check=="w"
#             print("warn")
#             @error("$w_msg $t : \n $data \n")
#         end


#         check, data, e_msg, w_msg = test_vals([zmort, dZdt, d_gain_total])
#         if check=="e"
#             @error("$e_msg $t : \n $data \n")
#         elseif check=="w"
#             print("warn")
#             @error("$w_msg $t : \n $data \n")
#         end


#         check, data, e_msg, w_msg = test_vals([dDdt])
#         if check=="e"
#             @error("$e_msg $t : \n $data \n")
#         elseif check=="w"
#             print("warn")
#             @error("$w_msg $t : \n $data \n")
#         end

import REPL
using REPL.TerminalMenus
using Printf, NCDatasets
using Dates
using SparseArrays



function message(v::String, nd::Int64=0, nb::Int64=0, np::Int64=0, nz::Int64=0, fsaven::String="")

    m = Dict(
        "START" => "\n -------------------------------- STARTING PROGRAM ----------------------------------- \n",
        "DEF" => """Default values: \n tt = 2 \n nrec = 100 \n nn = 1 \n np = 1 \n nz = 2 \n nb = 6 \n nd = 3 
                    pulse freq = 0 \n rsource = 0.1 \n rsink = 1.0 \n y_i = conts. 0.3 \n SW = const. 1.0 \n vmax_i = ordered \n """,
        "DF1" => ["Use defaults", "Select custom params"],
        "DF2" => "Proceed with defaults or select custom params?",
        "T" => "\nEnter simulation run time (tt - number of days): ",
        "REC" => "Enter number of timepoints to record (nrec): ",
        "OM" => "Enter number of organic matter pools (nd): ",
        "BA" => "Enter number of bacteria populations (nb): ",
        "NPZ" => """\n 
        -----------------------------------
        Currently running with default of: 
        1 phyto (np) \n  2 zoo (nz) \n  1 inorganic matter pool (nn) 
        -----------------------------------""",
        "CM" => "\n >> Building consumption matrix with $nd OM pools and $nb bac populations << ",
        "GM" => "\n >> Building grazing matrix with $np phy, $nb bac and $nz zoo populations << ",
        "Y1" => ["Equal for all bacteria", "Randomised"],
        "Y2" => "\n Select yield rates for bacteria populations (y_i):",
        "SW1" => ["Equal distribution", "Randomised"],
        "SW2" => "\n Select supply weight for OM pools (SW):",
        "SUB" => "\n SETTING SUBSTRATE TRAITS \n -------------------------- ",
        "UP1" => ["Ordered assignment", "Randomly selected along log range"],
        "UP2" => "Select max uptake rate (vmax_i, 1/d):",
        "S1"   => "\nEnter OM supply rate (mmol C/m3, decimal): ",
        "S2"   => "\nEnter OM sinking rate (1/day, decimal): ",
        "MIC" => "\n SETTING MICROBIAL TRAITS \n --------------------------",
        "GA1" => ["Tradeoff", "Constant affinity"],
        "GA2" => "Apply bacterial growth rate-affinity tradeoff?",
        "ENV" => "\n SETTING NUTRIENT SUPPLY \n -------------------------- ",
        "PU2" => "Simulate winter or summer conditions?",
        "PU1" => ["Winter", "Summer"],
        "SV" => "Saving to: $fsaven",

    )

    return m["$v"]

end


function set_savefiles(launch_time, years, nn, np, nz, nb, nd, nv, lysis=0)

    fsave = "results/outfiles/"

    if lysis == 0
        fsaven = string(fsave, Dates.format(launch_time, "yymmdd_HH:MM"), "_$(years)y_$(nn)N$(np)P$(nz)Z$(nb)B$(nd)D$(nv)V.nc")
    else
        fsaven = string(fsave, Dates.format(launch_time, "yymmdd_HH:MM"), "_$(years)y_$(nn)N$(np)P$(nz)Z$(nb)B$(nd)D$(nv)V_L.nc")
    end

    return fsaven

end


function get_defaults()

    tt = 5
    nrec = 100 
    nn = 1
    np = 1
    nz = 2
    nb = 6
    nd = 3
    y_i = ones(nd)*0.3 
    SW = ones(nd) 
    vmax_i = ordered_vmax(nd)
    pulse = 1
    rsource = 0.1
    rsink = 1.0

    return tt, nrec, nd, nb, np, nz, nn, y_i, SW, vmax_i, rsource, rsink, pulse

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


function update_tracking_arrs(track_n, track_p, track_z, track_b, track_d, track_v, track_time, ntemp, ptemp, ztemp, btemp, dtemp, vtemp, t, trec, prms)

    j = Int(t÷trec + 1)
    t_id = t.*prms.dt
    track_p[:,j] .= ptemp
    track_b[:,j] .= btemp 
    track_z[:,j] .= ztemp 
    track_n[:,j] .= ntemp 
    track_d[:,j] .= dtemp
    track_v[:,j] .= vtemp
    track_time[j] = t_id 

    @printf("Day %7.1f out of %5.0f = %4.0f%% done at %s \n", t_id, prms.days, t_id/prms.days*100, now())

    return track_n, track_p, track_z, track_b, track_d, track_v, track_time

end


function savetoNC(fsaven, p, b, z, n, d, v, timet, tst, tfn, prms, pulse)

    outdir = "/home/lee/Dropbox/Development/NPZBDV_0D/"
    path = joinpath(outdir, fsaven) 
    println("Saving to: ", path)

    # if pulse == 0
    #     season = "winter"
    # else
    #     season = "summer"
    # end

    f = NCDataset(path, "c") #c for create

    # define the dim of p, b, z, n, d
    defDim(f,"np", prms.np)
    defDim(f,"nb", prms.nb)
    defDim(f,"nz", prms.nz)
    defDim(f,"nn", prms.nn)
    defDim(f,"nd", prms.nd)
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
     w.attrib["units"] = "mmol/m3 C biomass"
 
     w = defVar(f,"b",Float64,("nb","nrec"))
     w[:,:] = b
     w.attrib["units"] = "mmol/m3 C biomass"
 
     w = defVar(f,"z",Float64,("nz","nrec"))
     w[:,:] = z
     w.attrib["units"] = "mmol/m3 C biomass"
     
     w = defVar(f,"n",Float64,("nn","nrec"))
     w[:,:] = n
     w.attrib["units"] = "mmol/m3 C OM"
 
     w = defVar(f,"d",Float64,("nd","nrec"))
     w[:,:] = d
     w.attrib["units"] = "mmol/m3 C OM"

     w = defVar(f,"v",Float64,("nv","nrec"))
     w[:,:] = v
     w.attrib["units"] = "mmol/m3 C OM"
     
    #  w = defVar(f,"uptake",Float64,("nd","nb"))
    #  w[:,:] = uptake
    #  w.attrib["units"] = "mmol/m3 C per d; uptake matrix"
     
     # w = defVar(f,"pIC",Float64,("np",))
     # w[:,:] = prms.pIC
     # w.attrib["units"] = "mmol/m3 C biomass"
 
     # w = defVar(f,"bIC",Float64,("nb",))
     # w[:,:] = prms.bIC
     # w.attrib["units"] = "mmol/m3 C biomass"
 
     # w = defVar(f,"zIC",Float64,("nz",))
     # w[:,:] = prms.zIC
     # w.attrib["units"] = "mmol/m3 C biomass"
 
     # w = defVar(f,"nIC",Float64,("nn",))
     # w[:,:] = prms.nIC
     # w.attrib["units"] = "mmol/m3 C OM"
 
     # w = defVar(f,"dIC",Float64,("nd",))
     # w[:,:] = prms.dIC
     # w.attrib["units"] = "mmol/m3 C OM"
     
     w = defVar(f, "timet", Float64, ("nrec",))
     w[:] = timet
     w.attrib["units"] = "days"
 
    #  w = defVar(f, "K_n", Float64, ("np",))
    #  w[:] = prms.K_n
    #  w.attrib["units"] = "m3/mmol; half-sat rate of p"
 
     w = defVar(f,"CM",Float64,("nd","nb"))
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

function define_dims(ds, prms, nrec1)

    vars = Dict("nb" => prms.nb, "nd" => prms.nd, "nrec" => nrec1)

    for (k, v) in x
        defDim(ds, k, v)
    end

    return ds

end

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
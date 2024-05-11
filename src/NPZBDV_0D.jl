module NPZBDV_0D
    
    import REPL
    using REPL.TerminalMenus
    using Logging, LoggingExtras, Dates
    using SparseArrays, Distributions

    # using Distributed
    # addprocs(14, exeflags = "--project=$(Base.active_project())")
    # println("Number of cores: ", nprocs())
    # println("Number of workers: ", nworkers())

    include("model.jl")
    include("params.jl")
    include("utils.jl")
    include("consumption_matrix.jl")
    include("grazing_matrix.jl")
    include("traits.jl")
    include("integrate.jl")
    include("plot.jl")


    #------------------------------------------------------------------------------------------------------------#
    #   SELECT LYSIS & GRAZING TYPE 
    #------------------------------------------------------------------------------------------------------------#
        lysis = request(message("LY1"), RadioMenu(message("LY2")))
        graze = request(message("GZ1"), RadioMenu(message("GZ2")))
        carbon = request(message("CB1"), RadioMenu(message("CB2")))
        pulse = 0

        rsource = 0.5
        rsink = 0.4
        # rsource = 0.0
        # rsink = 0.0
        CNr = 6

    #------------------------------------------------------------------------------------------------------------#
    #   SELECT SIMULATION TIME  
    #------------------------------------------------------------------------------------------------------------#
        println(message("TM"))
        user_input = readline()
        years = parse(Int64, user_input)

        days = years * 366
        nrec = days * 20
        dt = 0.01
        nt = Int(days/dt)

    #------------------------------------------------------------------------------------------------------------#
    #   LOAD MODEL
    #------------------------------------------------------------------------------------------------------------#
        files = readdir("src/prescribed")
        f = request("\nSelect model to run:", RadioMenu(files))
        include("prescribed/$(files[f])")
    
    #------------------------------------------------------------------------------------------------------------#
    #   RUN MODEL
    #------------------------------------------------------------------------------------------------------------#
        N, C, P, Z, B, Dn, Dc, V, track_time, fsaven = run_NPZBDV(params, lysis)
        log_params(params, lysis)
        plot_results(fsaven, lysis, graze, carbon, pulse)

end

export NPZBDV_0D


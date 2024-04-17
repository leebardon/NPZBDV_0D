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
    #   SELECT LYSIS TYPE 
    #------------------------------------------------------------------------------------------------------------#
        lysis = request(message("LY1"), RadioMenu(message("LY2")))
        pulse = 0

    #------------------------------------------------------------------------------------------------------------#
    #   SELECT SIMULATION TIME  
    #------------------------------------------------------------------------------------------------------------#
        simulation_time = request(message("TM1"), RadioMenu(message("TM2")))
        if simulation_time == 1
            years = 2
            days = 732
            nrec = 14640
        elseif simulation_time == 2
            years = 10
            days = 3660
            nrec = 73200
        elseif simulation_time == 3
            years = 50
            days = 18300
            nrec = 366000
        else 
            years = 100
            days = 36600
            nrec = 732000
        end

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
        N, P, Z, B, D, V, track_time, fsaven = run_NPZBDV(params, lysis)
        log_params(params, lysis)

        outdir = "/home/lee/Dropbox/Development/NPZBDV_0D/"
        plot_results(outdir, fsaven, lysis, pulse)

end

export NPZBDV_0D


    include("../model.jl")
    include("../params.jl")
    include("../utils.jl")
    include("../consumption_matrix.jl")
    include("../grazing_matrix.jl")
    include("../traits.jl")
    include("../integrate.jl")
    include("../plot.jl")


    #------------------------------------------------------------------------------------------------------------#
    #   CORE PARAMS 
    #------------------------------------------------------------------------------------------------------------#
        nn = 1
        nc = 1
        np = 6
        nz = 2
        nb = 2
        ndn = 2
        ndc = 2
        nv = 2  

        fsaven = set_savefiles(now(), years, nn, nc, np, nz, nb, ndn, ndc, nv, lysis)
        logger = set_logger(now())

    # -----------------------------------------------------------------------------------------------------------#
    #                                       PHYTOPLANKTON PARAMS 
    #------------------------------------------------------------------------------------------------------------
        CMp = [ 1 1 1 1 1 1 ]

        Fg_p = [0.1, 0.25, 0.5, 0.68, 0.79, 0.91]      # fraction of proteome optimized to growth
        Fa_p = 1. .- Fg_p                              # fraction optimized to substrate affintiy
    
        vmax_i = [0.5, 1.0, 2.0, 3.0, 4.0, 6.0]        # max growth rates (per day)
        Kp_i = vmax_i./10                              # half saturation of P_i
    
        vmax_ij = set_vmax_ij(nn, np, vmax_i, Fg_p)    # growth rate of P_i on N
        Kp_ij = set_Kp_ij(nn, np, Fa_p, CMp, vmax_ij)  # half saturation of P_i on N
    
    
    # -----------------------------------------------------------------------------------------------------------#
    #                                    HETEROTROPHIC BACTERIA PARAMS
    #------------------------------------------------------------------------------------------------------------#
        CM = [1 1 
              0 0 ] 
        
        y_i = ones(ndn)*0.3
        y_ij = broadcast(*, y_i, CM) 
    
        Fg_b = [1.0, 1.0]
        Fa_b = 1. .- Fg_b      
    
        umax_i = [1.0, 1.0]      # Fg_b and umax_i are dummy vals - already provided from Emily's previous work (trade-off applied)
        Km_i = umax_i./10 
    
        umax_ij =  [32.0  4.0
                    0.0  0.0 ] 
    
        Km_ij = [1.54  0.124
                 0.0  0.0 ]
    
    # -----------------------------------------------------------------------------------------------------------#
    #                ZOOPLANKTON PARAMS 
    #------------------------------------------------------------------------------------------------------------#
        if graze == 1
            # GrM - first 6 cols are phyto, next 1 dom consuming bacteria, rows are zoo
            GrM = [ 1.0  1.0  1.0  1.0  1.0  1.0  0.0  0.0
                    0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0 
                ] ;
        
            g_max = ones(nz)
            K_g = ones(nz)*1.0
            γ = ones(nz)*0.3
            zIC = ones(nz)*0.01
        else
            GrM = fill(0, (nz, (np + nb)))
            g_max = zeros(nz)
            K_g = zeros(nz)
            γ = zeros(nz)
            zIC = zeros(nz) 
        end


    # -----------------------------------------------------------------------------------------------------------#
    #                                     MORTALITY RATES (mmol/m3/day)
    #------------------------------------------------------------------------------------------------------------#
        m_lp = ones(np) * 0.1
        m_qp = ones(np) * 0.1  # (.1 if explicit grazers, if not, 1)
    
        m_lb = ones(nb) * 0.01
        m_qb = ones(nb) * 0.1
    
        m_lz = ones(nz) * 0.01
        m_qz = ones(nz) * 1.0

    # -----------------------------------------------------------------------------------------------------------#
    #   ORGANIC MATTER
    #------------------------------------------------------------------------------------------------------------#
        if carbon == 1
            cIC = ones(nc)*20.0
            dcIC = ones(ndc)*0.1
        else
            cIC = ones(nc)*0.0
            dcIC = ones(ndc)*0.0
        end
    
        # distribution of OM from mortality and lysis to detritus pools
        # viral decay amost entirely contributes to labile pool, lysis weighted to labile but includes both (walls & innards)
        om_dist_mort = [0.6, 0.4]   
        om_dist_lys = [0.8, 0.2]   
        om_dist_vde = [1.0, 0.0]  
    

    # -----------------------------------------------------------------------------------------------------------#
    #                                       VIRUS PARAMS 
    #------------------------------------------------------------------------------------------------------------#
        if lysis == 1
            # virus for each B
            VM = [ 1  0
                   0  1 ]
            
            vly = 8e-15             # m3/virus/dat (from weitz et al 2015)
            vbs = 25                # burst size 
            vde = 0.17              # decay rate (day-1)  0.17
            vIC = ones(nv)*0.1
            # how much nitrogen per virus? get a virus quota per nitrogen i.e. mmol N / virus then convert back to mmol/day
        else 
            VM = fill(0, (ndn, nb))
            vly = 0.0
            vbs = 0                
            vde = 0.0   
            vIC = zeros(nv)        
        end

    # -----------------------------------------------------------------------------------------------------------#
    #   INITIAL CONDITIONS
    #------------------------------------------------------------------------------------------------------------#
        nIC = ones(nn)*5.0
        pIC = ones(np)*0.1
        bIC = ones(nb)*0.01
        dnIC = ones(ndn)*0.1

        dnIC[2] = 0.0

    # -----------------------------------------------------------------------------------------------------------#
    #   INITIATE PARAMS
    #------------------------------------------------------------------------------------------------------------#
    params = Prms(
                years, days, nrec, dt, nt, 
                nn, nc, np, nz, nb, ndn, ndc, nv, CNr, nIC, cIC, pIC, zIC, bIC, dnIC, dcIC, vIC,
                vmax_i, vmax_ij, Kp_i, Kp_ij, m_lp, m_qp, CMp, Fg_p,
                umax_i, umax_ij, Km_i, Km_ij, y_ij, m_lb, m_qb, CM, Fg_b,
                g_max, K_g, γ, m_lz, m_qz, GrM, vly, vbs, vde, VM, 
                fsaven, om_dist_mort, om_dist_lys, om_dist_vde, rsource, rsink, pulse, graze, lysis, carbon
            )


    #   RUN MODEL
    #------------------------------------------------------------------------------------------------------------#
    # @info("Model Params: \n $params \n")
    # N, P, Z, B, D, V, track_time, fsaven = run_NPZBDV(params, lysis)
    # log_params(params, lysis)

    # outdir = "/home/lee/Dropbox/Development/NPZBDV_0D/"
    # plot_results(outdir, fsaven, lysis, pulse)
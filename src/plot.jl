using NCDatasets
using Plots, ColorSchemes
using DataFrames

default(show = true)

function plot_results(outdir, outfile, lysis, pulse=0)

    ds = NCDataset(outdir*outfile)
    file = replace(outfile, "out_0D_" => "", ".nc" => "", "results/outfiles/" => "")
    pulse == 1 ? season = "win" : season = "sum"

    n, c, p, z, b, d, v, timet = get_plot_variables(ds, ["n", "c", "p", "z", "b", "d", "v", "timet"])

    f1 = plot_total(n, c, p, z, b, d, v, timet)
    pb, pp = lysing(p, b, lysis)
    # f2 = plot_all_pzbv(n, p, z, b, d, v, timet)
    f3 = plot_individual(n, c, p, z, b, d, v, timet)
    f4 = plot_total_pzbv(p, z, b, v, timet)

    savefig(f1,"results/plots/$(file)_1.pdf")
    # savefig(f2,"results/plots/$(file)_2.pdf")
    savefig(pb,"results/plots/lysing/$(file)_B.pdf")
    savefig(pp,"results/plots/lysing/$(file)_P.pdf")
    savefig(f3,"results/plots/$(file)_3.pdf")
    savefig(f4,"results/plots/$(file)_4.pdf")

end

function lysing(p, b, lysis)

    B = b[:, end];
    P = p[:, end];
    np, nb = get_size([p,b])
    xb = collect(1:nb);
    xp = collect(1:np);

    lysis == 1 ? tit="Explicit Lysing" : tit="Implicit Lysing"
    pb = plot(xb, B, title=tit, xlabel="B #", ylabel = "Conc.", label="B", grid=false, lw=5);
    pp = plot(xp, P, title=tit, xlabel="P #", ylabel = "Conc.", label="P", grid=false, lw=5, lc="seagreen"); 

    return pb, pp

end


function plot_total(n, c, p, z, b, d, v, timet)

    all_n, all_c, all_p, all_z, all_b, all_d, all_v = sum_subtypes([n, c, p, z, b, d, v])

    plt = plot(timet[:], [all_p[:], all_b[:], all_z[:], all_v[:]], lw=4, lc=[:limegreen :skyblue3 :coral3 :black], label=["P" "B" "Z" "V"], grid=false, alpha=0.7)
    plot!(timet[:], [all_n[:], all_c[:], all_d[:]], lw=1.5, linecolor=[:darkgrey :purple :olive], ls=:dashdot, label=["N" "DIC" "DOM"], alpha=1.0)

    title!("Total NCPZBDV over time")
    xlabel!("Time (days)")
    ylabel!("Concentration")

    return plt

end


function plot_all_pzbv(n, p, z, b, d, v, timet)

    nn, np, nz, nb, nd, nv = get_size([n,p,z,b,d,v])

    plt = plot(timet[:], p[1, :], lw=2, linecolor="limegreen", grid=false, label="P")

    if np > 1
        cp1 = Plots.palette(:greens, np)
        for i in 2:np
            plot!(timet[:], p[i, :], lw=3.5, palette=cp1, grid=false, label=" ", alpha=0.7)
        end
    end 

    cp2 = Plots.palette(:dense, nb)
    for j in 1:nb
        j == 1 ? lab="B" : lab=" "
        plot!(timet[:], b[j, :], lw=2.5, palette=cp2, grid=false, label=lab, alpha=0.8)
    end

    cp3 = Plots.palette(:solar, nz)
    for k in 1:nz
        k == 1 ? lab="Z" : lab=" "
        plot!(timet[:], z[k, :], lw=1.5, color="black", grid=false, label=lab)
    end

    cp4 = Plots.palette(:reds, nv)
    for x in 1:nv
        x == 1 ? lab="V" : lab=" "
        plot!(timet[:], v[x, :], lw=2.0, palette=cp4, grid=false, label=lab, ls=:dot)
    end

    # cp4 = Plots.palette(:turbid, nd)
    # for l in 1:nd
    #     l == 1 ? lab="D" : lab=" "
    #     plot!(timet[:], d[l, :], lw=1.5, ls=:dash, palette=cp3, grid=false, label=lab)
    # end

    title!("PZBV over time")
    xlabel!("Time (days)")
    ylabel!("Concentration (mmol/m^3)")

    return plt

end

function plot_individual(n, c, p, z, b, d, v, timet)

    all_n, all_c, all_p, all_z, all_b, all_d, all_v = sum_subtypes([n, c, p, z, b, d, v])

    plt = plot(timet[:], [all_p[:], all_b[:], all_z[:], all_v[:], all_d[:], all_n[:], all_c[:]]; 
    layout = 7, 
    linewidth = 3.5, 
    alpha=0.5,
    legend = false, 
    title = ["Total P" "Total B" "Total Z" "Total V" "Total DOM" "Total N" "Total DIC"]
    )

    return plt

end


function plot_total_pzbv(p, z, b, v, timet)

    all_p, all_z, all_b, all_v = sum_subtypes([p, z, b, v])

    plt = plot(timet[:], [all_p[:], all_b[:], all_z[:], all_v[:]], lw=5, lc=[:limegreen :skyblue3 :coral3 :black], label=["P" "B" "Z" "V"], 
    grid=false, alpha=0.7)
    # plot!(timet[:], all_d[:], lw=2, linecolor=[:darkgrey :olive], ls=:dash, label=["D"])

    title!("Total PZBV over time")
    xlabel!("Time (days)")
    ylabel!("Concentration")

    return plt
end


#---------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------


function sum_subtypes(all_states)

    out = Vector{Any}()

    for s in all_states
        append!(out, [sum(s[:,:], dims=1)])
    end

    return out
end


function get_plot_variables(ds, vars)

    out = Vector{Any}()

    for v in vars
        append!(out, [ds["$v"]])
    end

    return out

end

function get_size(arr)

    out = Vector{Int}()
    
    for a in arr
        append!(out, size(a,1))
    end

    return out

end


# plot_results(outdir, outfile, 1)
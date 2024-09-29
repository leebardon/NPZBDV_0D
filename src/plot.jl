using NCDatasets
using Plots, LaTeXStrings, ColorSchemes
using DataFrames

# default(show = true)
include("utils.jl")
outdir = "/home/lee/Dropbox/Development/NPZBDV_0D/"

function plot_results(outfile, lysis, graze, carbon, pulse=0)

    ds = NCDataset(outdir*outfile)
    file = replace(outfile, "out_0D_" => "", ".nc" => "", "results/outfiles/" => "")

    n, c, p, z, b, dn, dc, v, timet = get_plot_variables(ds, ["n", "c", "p", "z", "b", "dn", "dc", "v", "timet"])

    f1 = plot_substrates(n, c, p, z, b, dn, dc, v, timet)
    f2 = plot_heterotrophs(p, z, b, v, timet)
    f3 = plot_bio(p, z, b, v, timet)
    f4 = plot_individual(n, c, p, z, b, dn, dc, v, timet)

    l = @layout [
        a{0.33w} b{0.33w} c{0.33w} 
                 d{0.6h}  
    ]

    f = plot(f1, f2, f3, f4,
    fg_legend = :transparent,
    layout = l,
    )
    savefig(f,"results/plots/$(file).png")

    # if carbon == 1
    #     f5 = plot_OM(dn, dc, timet)
    #     savefig(f5,"results/plots/$(file)_OM.png")
    # end

    # pb, pp = lysing(p, b, lysis)
    # savefig(pb,"results/plots/lysing/$(file)_B.pdf")
    # savefig(pp,"results/plots/lysing/$(file)_P.pdf")

end

function plot_structure(filename, ds1, ds2, ds3=0, ds4=0)

    B1 = get_endpoints(["b"], ds1)
    B2 = get_endpoints(["b"], ds2)
    al=0.6

    # scatter(B1, grid=false, markersize=10, alpha=al, label="Implicit", show=true, fg_legend=:transparent,
    #         legendfontsize=10, xlabel="Microheterotroph", ylabel=L"mmol N/m^3");
    # scatter!(B2, grid=false, markersize=10, label="Lysis", alpha=al);

    # if ds3 != 0
    #     B3 = get_endpoints(["b"], ds3)
    #     scatter!(B3, grid=false, markersize=10, label="Grazers", alpha=al);
    # end

    plot(B1, grid=false, lw=8, alpha=al, label="Implicit", show=true, fg_legend=:transparent,
            legendfontsize=10, xlabel="Microheterotroph", ylabel=L"mmol N/m^3");
    plot!(B2, grid=false, lw=10, label="Lysis", alpha=al);

    if ds3 != 0
        B3 = get_endpoints(["b"], ds3)
        plot!(B3, grid=false, lw=8, label="Grazers", alpha=al);
    end

    if ds4 != 0
        B4 = get_endpoints(["b"], ds4)
        plot!(B4, grid=false, lw=8, label="L & G", alpha=al);
    end


    # fig2 = plot(timet[:], b[1, :], lw=l, linecolor="seagreen", label=false, alpha=alp)
    # if nb > 1
    #     for i in 2:nb
    #         plot!(timet[:], b[i, :], lw=l, palette=:batlow10, label=false, alpha=alp)
    #     end
    # end 

    title!("Microhet. Community Structure")

    savefig(filename)

end 

function plot_substrates(n, c, p, z, b, dn, dc, v, timet)

    all_n, all_c, all_p, all_z, all_b, all_dn, all_dc, all_v = sum_subtypes([n, c, p, z, b, dn, dc, v])

    alp=0.8
    start_ts=100
    
    fig1 = plot(timet[:], [all_n[:], all_dn[:]], lw=3, linecolor=[:orange :darkred], ls=:dot, label=[" DIN" " DON"])

    if sum(all_c) > 0
        plot!(timet[:], [all_c[:], all_dc[:]], lw=3, linecolor=[:grey82 :grey30], ls=:dot, label=[" DIC" " DOC"])
    end

    title!("Substrate Dynamics")
    # xlabel!("Time (days)")
    ylabel!(L" mmol ~N/m^3")

    return fig1

end

function plot_heterotrophs(p, z, b, v, timet)

    np, nz, nb, nv = get_size([p,z,b,v])
    alp=0.7
    l=4

    fig2 = plot(timet[:], b[1, :], lw=l, linecolor="seagreen", label=false, alpha=alp)
    if nb > 1
        for i in 2:nb
            plot!(timet[:], b[i, :], lw=l, palette=:batlow10, label=false, alpha=alp)
        end
    end 

    title!("Microhet. Dynamics")

    return fig2

end

function plot_bio(p, z, b, v, timet)

    np, nz, nb, nv = get_size([p,z,b,v])
    alp=0.7
    l=4

    fig3 = plot(timet[:], p[1, :], lw=l, linecolor="seagreen", label=false, alpha=alp)
    if np > 1
        for i in 2:np
            plot!(timet[:], p[i, :], lw=l, color="seagreen", label=false, alpha=alp)
        end
    end 

    if nz > 0
        for k in 1:nz
            k == 1 ? lab="Z" : lab=" "
            plot!(timet[:], z[k, :], lw=l, color="black", label=false, alpha=alp)
        end
    end

    if nv > 0
        for x in 1:nv
            x == 1 ? lab="V" : lab=" "
            plot!(timet[:], v[x, :], lw=l, color="purple", label=false, alpha=alp)
        end
    end

    title!("PZV dynamics")

    return fig3

end


function plot_individual(n, c, p, z, b, dn, dc, v, timet)

    all_n, all_c, all_p, all_z, all_b, all_dn, all_dc, all_v = sum_subtypes([n, c, p, z, b, dn, dc, v])

    fig = Array{Plots.Plot, 1}(undef, 8);
    lsize = 4
    lg=false

    fig[1] = plot(timet[:], all_n[:], lw=lsize, legend=lg, ylabel=L" mmol ~N/m^3", title="Total N") 
    fig[2] = plot(timet[:], all_dn[:], lw=lsize, legend=lg, ylabel=" ", title="Total DON") 
    fig[3] = plot(timet[:], all_p[:], lw=lsize, legend=lg, ylabel=" ", title="Total Phy") 
    fig[4] = plot(timet[:], all_b[:], lw=lsize, legend=lg, ylabel=" ", title="Total Bac")

    if sum(all_z[:]) > 0
        fig[5] = plot(timet[:], all_z[:], lw=lsize, legend=lg, ylabel=L" mmol ~N/m^3", title="Total Zoo") 
    else 
        fig[5] = plot(timet[:], all_z[:], lw=0, legend=lg, ylabel=L" mmol ~N/m^3", ylim=(0,5), grid=false, title="Z (none)") 
    end

    if sum(all_v[:]) > 0
        fig[6] = plot(timet[:], all_v[:], lw=lsize, legend=lg, ylabel=" ", title="Total Vir") 
    else 
        fig[6] = plot(timet[:], all_v[:], lw=0, legend=lg, ylabel=" ", ylim=(0,5), grid=false, title="V (none)") 
    end

    if sum(all_c[:]) > 0
        fig[7] = plot(timet[:], all_c[:], lw=lsize, legend=lg, ylabel=" ", title="Total C") 
        fig[8] = plot(timet[:], all_dc[:], lw=lsize, legend=lg, ylabel=" ", title="Total DOC")
    else 
        fig[7] = plot(timet[:], all_c[:], lw=0, legend=lg, ylabel=" ", ylim=(0,5), grid=false, title="C (none)") 
        fig[8] = plot(timet[:], all_dc[:], lw=0, legend=lg, ylabel=" ", ylim=(0,5), grid=false, title="DOC (none)")
    end

    f = plot(fig..., 
    fg_legend = :transparent,
    layout = (2,4),
    size=(800,800),
    )

    return f
end


function plot_OM(dn, dc, timet)

    ndn, ndc = get_size([dn, dc])

    plt = plot(timet[:], dn[1, :], lw=2, linecolor="limegreen", grid=false, label="DON")

    cp1 = Plots.palette(:thermal, dn)
    for j in 2:ndn
        plot!(timet[:], dn[j, :], lw=4, palette=cp1, grid=false, label=" ")
    end

    cp2 = Plots.palette(:turbid, dc)
    for l in 1:ndc
        l == 1 ? lab="DOC" : lab=" "
        plot!(timet[:], dc[l, :], lw=4, ls=:dot, palette=cp2, grid=false, label=lab)
    end

    title!("DON & DOC")
    xlabel!("Time (days)")
    ylabel!(L" mmol ~N/m^3")

    return plt

end


function lysing(n, p, z, b, d, v, lysis)

    B, P = b[:, end], p[:, end];
    np, nb = get_size([p,b])
    xb, xp = collect(1:nb), collect(1:np);
    lysis == 1 ? tit="Explicit Lysing" : tit="Implicit Lysing"

    fig = Array{Plots.Plot, 1}(undef, 2);
    fig[1] = plot(xb, B, xlabel="B #", ylabel=L" mmol ~N/m^3", label="Bac", grid=false, lw=5);
    fig[2] = plot(xp, P, xlabel="P #", yformatter=Returns(""), label="Phy", grid=false, lw=5, lc="seagreen"); 

    f = plot(fig..., 
    fg_legend = :transparent,
    layout = (1,2),
    size=(400,300),
    title=tit,
    )

    return f

end

function lysis_compare(with_ds, without_ds)

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



# fsaven1 = "results/outfiles/240604_15:36_10y_1N1C6P2Z5B10D5V.nc";
# fsaven2="results/outfiles/240604_15:34_10y_1N1C6P2Z5B10D5V_L.nc";
# ds1=NCDataset(fsaven1);
# ds2=NCDataset(fsaven2);
# plot_structure(ds1, ds2)

fsaven1 = "results/outfiles/240902_13:28_50y_1N1C6P2Z5B10D5V.nc";
fsaven2 = "results/outfiles/240902_14:24_20y_1N1C6P2Z5B10D5V_L.nc";
fsaven3 = "results/outfiles/240902_13:35_50y_1N1C6P2Z5B10D5V.nc";
fsaven4 = "results/outfiles/240902_14:33_50y_1N1C6P2Z5B10D5V_L.nc";
ds1 = NCDataset(fsaven1);
ds2 = NCDataset(fsaven2);
ds3 = NCDataset(fsaven3);
ds4 = NCDataset(fsaven4);
plot_structure("new_5B_quad2.png", ds1, ds2, ds3, ds4)

# fsaven = "results/outfiles/240902_12:41_50y_1N1C6P2Z5B10D5V_L.nc";
# lysis, graze, carbon = 1, 2, 2
# plot_results(fsaven, lysis, graze, carbon);
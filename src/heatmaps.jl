# using NCDatasets, DataFrames
# using CairoMakie
# using LinearAlgebra, Statistics

# include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/utils.jl")
# include("/home/lee/Dropbox/Development/NPZBD_1D/src/utils/save_utils.jl")

function get_final_year(ds, vars)
    # where year is 12 * 30 days and 20 timesteps recorded per day
    ts = 20 * 366

    final = Vector{Any}()

    for v in vars
        append!(final, [ds[v][:, end-(ts-1):end]])
    end

    return final

end;


function set_zmax(state_var, num_state_var)

    zmax, z97 = 0, 0

    for s in range(1, num_state_var)
        state_var_max = maximum(state_var[:, s, :])
        state_var_max > zmax ? zmax = state_var_max : nothing

        s_arr = vec(state_var[:, s, :])
        state_var_97 = quantile(s_arr, 0.97)
        state_var_97 > z97 ? z97 = state_var_97 : nothing
    end

    return zmax, z97

end


function plot_bmass_heatmaps(fsaven)

    ds = NCDataset(fsaven)
    filename = replace(fsaven, "_L"=> "", ".nc" => "", "/home/lee/Dropbox/Development/NPZBDV_0D/" => "", "results/outfiles/" => "")

    p, b = get_final_year(ds, ["p", "b"])

    # bmass_heatmaps(p, filename, "P")
    bmass_heatmaps(b, filename, "B")

end


function bmass_heatmaps(state_var, filename, varname)

    lab = sum(state_var[[1, 11, 21, 2, 12, 22, 3, 13, 23, 4, 14, 24], :], dims=1)
    semi_lab = sum(state_var[[5, 15, 25, 6, 16, 26], :], dims=1)
    refract = sum(state_var[[7, 17, 27, 8, 18, 28, 9, 19, 29, 10, 20, 30], :], dims=1)

    lfs=9

    # joint_limits = (0, zmax)
    # row, col = 1, 1
    # fig = Figure(size=res)
    # for i in range(1, 3)
    #     data = state_var[1:max_d, i]
    #     survivors = set_extinct_to_zero(data)
    #     daily_data = survivors[:, 1:20:end]
    #     days = collect(1:size(daily_data, 2))
    #     z = get_hmap_z_axis(d, days, daily_data)

    #     heatmap(fig[row, col], days, reverse(d), reverse(z), xrotation=45,  colorrange=joint_limits)
    #     col += 1
    #     if mod(i, 5) == 0 
    #         row += 1
    #         col -= 5
    #     end
    # end

    Colorbar(fig[:, end+1], colorrange=joint_limits, size=30, label=L"mmol/m^3", labelsize=20)

    save("results/heatmaps/$(varname)_$(filename).png", fig)
end


function npp_heatmaps(fsaven, npp_months, npp, mov_avs, daily_means)
       
    parent_folder = "results/plots/heatmaps/npp/"
    filename = replace(fsaven, ".nc" => "", "/home/lee/Dropbox/Development/NPZBD_1D/" => "", "results/outfiles/" => "")
    dir = check_subfolder_exists(filename, parent_folder)

    cols_hm = :tofino25
    cols_ln = ["orangered", "seagreen3", "purple4", "black", "skyblue1"]
    dnames = ["5m", "50m", "100m", "150m", "200m"]

    zc = get_zc(890)
    depth = -zc[1:30]
    ts = collect(1:length(mov_avs[1]))

    let
        fig = Figure(resolution=(600,1250))

        z = npp'
        limits = (minimum(npp), maximum(npp))
        ax1, hm = heatmap(fig[1:2, 1], ts, reverse(depth), reverse(z), clim=limits, colormap=(cols_hm))
        Colorbar(fig[1:2, 2], colorrange=limits, colormap=(cols_hm), size=20, label=L"mmol/m^3", labelsize=20)
        ax1.xlabel=L"Timesteps"
        ax1.ylabel=L"Depth~(m)"
        ax1.title="Net Primary Productivity"

        z2 = npp_months'
        days = collect(1:size(npp_months, 2))
        limits = (minimum(npp_months), maximum(npp_months))
        ax2, hm = heatmap(fig[3:4, 1], days, reverse(depth), reverse(z2), clim=limits, colormap=(cols_hm))
        Colorbar(fig[3:4, 2], colorrange=limits, colormap=(cols_hm), size=20, label=L"mmol/m^3", labelsize=20)
        ax2.xlabel=L"Months"
        ax2.xticks = 1:12
        ax2.ylabel=L"Depth~(m)"
        ax2.title="NPP Monthly Means"

        ax3 = fig[5, 1] = Axis(fig)
        for i in range(1, length(dnames))
            lines!(ts, mov_avs[i], alpha=0.7, linewidth=3, color=cols_ln[i], label=dnames[i])
             # td = collect(1:length(daily_means[1]))
            # scatterlines!(t, convert(Array{Float64}, daily_means[i]), markersize=3, color=cols_ln[i], strokewidth=1, strokecolor=:black, alpha=0.1, label=dnames[i])
            # lines!(td, convert(Array{Float64}, daily_means[i]), color=cols_ln[i], alpha=0.5, label=dnames[i], grid=false)
        end

        fig[5, 2] = Legend(fig, ax3, framevisible=false)
        ax3.xlabel=L"Timesteps"
        ax3.ylabel=L"mmol/m^3"
        ax3.title="Moving Averages"

        println("Saving fig to $(dir)/$(filename)_npp.png")
        save("$(dir)/$(filename)_npp.png", fig)
    end

end


fsaven = "results/outfiles/240428_14:42_10y_1N6P2Z30B10D30V_L.nc"
# fsaven = "results/outfiles/240428_14:46_10y_1N6P2Z30B10D30V_L.nc"

plot_bmass_heatmaps(fsaven)
# plot_bmass_heatmaps(fsaven, "B")


# plot_bmass_heatmaps(fsaven, "B")

# for var in ["P", "B"]
#     plot_bmass_heatmaps(fsaven, var, bloom)
# end
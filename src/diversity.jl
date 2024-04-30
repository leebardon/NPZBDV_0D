# b_imp_tot = sum(b_imp);
# b_imp_RA = b_imp ./ b_imp_tot;

# p3 = plot(nb, b_imp_RA, lw=5, alpha=0.7, grid=false, title="Relative Abundance after 10 years", 
# xlabel="Heterotrophic Bacteria", ylabel="RA", foreground_color_legend=nothing, label="Implicit" );
# plot!(nb, b_same_RA, lw=5, alpha=0.7, label="Explicit v1")
# plot!(nb, b_diff_RA, lw=5, alpha=0.7, label="Explicit v2")
# plot!(nb, b_exp_nz_RA, lw=5, alpha=0.7, label="Explicit v2 (no Z)")

# shan_imp = -sum(b_imp_RA .* log.(b_imp_RA))

# p4 = scatter([1], [shan_imp], markersize=8, grid=false, title="Shannon Diversity after 10 years", 
# xlabel="Group", ylabel="Shannon Diversity", label="Implicit", foreground_color_legend=nothing, xlims=(0, 4.5), 
# ylims=(3.1,3.4))
# scatter!([2], [shan_same], markersize=8, label="Explicit v1")
# scatter!([3], [shan_diff], markersize=8, label="Explicit v2")
# scatter!([4], [shan_exp_nz], markersize=8, label="Explicit v2 (no Z")


# lab=[1, 11, 21, 2, 12, 22, 3, 13, 23, 4, 14, 24];
# slab=[5, 15, 25, 6, 16, 26];
# sref=[7, 17, 27, 8, 18, 28, 9, 19, 29, 10, 20, 30];

# b_imp_lab_RA = sum(b_imp[lab]) / b_imp_tot
# b_imp_slab_RA = sum(b_imp[slab]) / b_imp_tot
# b_imp_sref_RA = sum(b_imp[sref]) / b_imp_tot
# b_same_lab_RA = sum(b_same[lab]) / b_same_tot
# b_same_slab_RA = sum(b_same[slab]) / b_same_tot
# b_same_sref_RA = sum(b_same[sref]) / b_same_tot
# b_diff_lab_RA = sum(b_diff[lab]) / b_diff_tot
# b_diff_slab_RA = sum(b_diff[slab]) / b_diff_tot
# b_diff_sref_RA = sum(b_diff[sref]) / b_diff_tot
# b_exp_nz_lab_RA = sum(b_exp_noz[lab]) / b_exp_nz_tot
# b_exp_nz_slab_RA = sum(b_exp_noz[slab]) / b_exp_nz_tot
# b_exp_nz_sref_RA = sum(b_exp_noz[sref]) / b_exp_nz_tot

# scatter(["Implicit"], [b_imp_lab_RA], xlims=(0, 4.5), xrotation=45, grid=false, 
# foreground_color_legend=nothing, markersize=8, color=:maroon, title="RA by Functional Type", 
# ylabel="RA", label="Lab.")
# scatter!(["Implicit"], [b_imp_slab_RA], color=:grey, label="Semi-Lab.", markersize=8)
# scatter!(["Implicit"], [b_imp_sref_RA], color=:darkcyan, label="Semi-Ref.", markersize=8)
# scatter!(["Exp. v1"], [b_same_lab_RA], color=:maroon, label="", markersize=8)
# scatter!(["Exp. v1"], [b_same_slab_RA], color=:grey, label="", markersize=8)
# scatter!(["Exp. v1"], [b_same_sref_RA], color=:darkcyan, label="", markersize=8)
# scatter!(["Exp. v2"], [b_diff_lab_RA], color=:maroon, label="", markersize=8)
# scatter!(["Exp. v2"], [b_diff_slab_RA], color=:grey, label="", markersize=8)
# scatter!(["Exp. v2"], [b_diff_sref_RA], color=:darkcyan, label="", markersize=8)
# scatter!(["Exp. no Z"], [b_exp_nz_lab_RA], color=:maroon, label="", markersize=8)
# scatter!(["Exp. no Z"], [b_exp_nz_slab_RA], color=:grey, label="", markersize=8)
# scatter!(["Exp. no Z"], [b_exp_nz_sref_RA], color=:darkcyan, label="", markersize=8)




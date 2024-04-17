using Statistics, StatsBase, Random


function build_grazing_matrix(nb, np, nz)

    GrM = collect(gmatrix(nb, np, nz))

    return GrM
end


function gmatrix(nb, np, nz)
 
    GM = sprand(Bool, nz, (nb+np), 0.6)
    GrM = check_for_empty_cols(GM)

    return GrM

end

function check_for_empty_cols(GM)

    empty_cols = findall(x -> x == 0, sum(GM, dims=1))
    s = size(empty_cols)

    if s[1] > 0
        for i in 1:s[1]
            new_col = sprand(Bool, nz, 1, 0.6)
            GM[:, empty_cols[i][2]] = new_col
        end
        GM = check_for_empty_cols(GM)
    end

    return GM

end



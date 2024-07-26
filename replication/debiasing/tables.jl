using CSV
using DataFrames
using Printf

function lifetime_method_order(l)
    if l == "opt"
        return 2
    elseif l == "pol"
        return 1
    end
end

data = CSV.read("./replication/debiasing/results.csv", DataFrame)
data = select!(data, sort(names(data)))
#data = select!(data, ["d", "n", "B", "J_estimator", "J_lifetime", "lifetime_method",
#"lifetime_multiplier", "lambda", "rmse", "bias", "sd",
#"bias_over_sd", "sd_hat",
#"sigma2_hat", "bias_theory", "sd_theory", "coverage", "average_width"])

data = sort!(data, [:d, :n, :B, :J_estimator, :J_lifetime, order(:lifetime_method, by=lifetime_method_order),
                    :lifetime_multiplier])

function make_table(df)
    d = df[1, "d"]
    n = df[1, "n"]
    B = df[1, "B"]
    tex = "\\begin{tabular}{|cc|cc|cccc|cc|cc|cc|}\n"
    tex *= "%\$d=$d\$, & \$n=$n\$, & \$B=$B\$&&&&&&&&&&\\\\\n"
    tex *= "\\hline\n"
    tex *= "\$J\$ & LS & LM & \$\\lambda\$ & RMSE & Bias & SD & Bias/SD & "
    tex *= "\$\\widehat{\\textrm{SD}}\$ & \$\\hat\\sigma^2\$ & ABias & ASD & CR & CIW \\\\\n"

    #display(df)
    for i in 1:nrow(df)
        row = df[i, :]

        if i > 1 && df[i, :J_lifetime] == df[i-1, :J_lifetime] &&
                     df[i, :J_estimator] == df[i-1, :J_estimator]
        else
            tex *= "\\hline\n"
        end

        if i > 1 && df[i, :J_estimator] == df[i-1, :J_estimator]
            tex *= ""
        else
            tex *= "$(df[i, :J_estimator])"
        end

        if i > 1 && df[i, :J_lifetime] == df[i-1, :J_lifetime] &&
                     df[i, :lifetime_method] == df[i-1, :lifetime_method]
            tex *= "&"
        else
            Jl = df[i, :J_lifetime]
            lm = df[i, :lifetime_method]
            if lm != "opt"
                hat = "\\hat"
            else
                hat = ""
            end
            tex *= "& \$$hat\\lambda_{$Jl}^{\\scriptsize{\\textrm{$lm}}}\$"
        end

        for col in [:lifetime_multiplier, :lambda, :rmse, :bias, :sd,
                    :bias_over_sd, :sd_hat, :sigma2_hat, :bias_theory, :sd_theory,
                    :coverage, :average_width]
            cell = df[i, col]
            if col == :coverage
                cell = 100 * cell
                cell = @sprintf "%.1f" cell
                cell = "$cell\\%"
            elseif col == :lifetime_multiplier
                cell = @sprintf "%.1f" cell
            elseif isa(cell, Float64)
                cell = @sprintf "%.4f" cell
            end
            tex *= "& $cell"
        end

        tex *= "\\\\\n"
    end
    tex *= "\\hline\n"
    tex *= "\\end{tabular}"
    write("./replication/debiasing/table_d$(d)_n$(n)_b$B.tex", tex)
end

for d in unique(data[!, "d"])
    println(d)
    data_d = filter(:d => ==(d), data)
    for n in unique(data_d[!, "n"])
        println(n)
        data_d_n = filter(:n => ==(n), data_d)
        make_table(data_d_n)
        display(data_d_n)
    end
end

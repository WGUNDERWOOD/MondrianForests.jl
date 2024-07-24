using CSV
using DataFrames

function lambda_method_order(l)
    if l == "optimal"
        return 1
    elseif l == "polynomial"
        return 2
    elseif l == "gcv"
        return 3
    end
end

function lambda_method_format(l)
    if l == "optimal"
        return "OPT"
    elseif l == "polynomial"
        return "POLY"
    elseif l == "gcv"
        return "GCV"
    end
end

data = CSV.read("./replication/debiasing/results.csv", DataFrame)
data = select!(data, sort(names(data)))
data = select!(data, ["d", "n", "B", "J_estimator", "J_lifetime", "lifetime_method",
                      "lifetime_multiplier", "lambda", "rmse", "bias", "sd_hat",
                      "sigma2_hat", "bias_theory", "sd_theory", "coverage", "average_width"])

# TODO work on this
#data = sort!(data, [:J_estimator, :lambda_target, order(:lambda_method, by=lambda_method_order),
                    #:J_lifetime, :d, :n, :B_estimator, :B_lifetime, :n_subsample,
                    #:lambda_multiplier])
#display(data)

function make_table(df)
    d = df[1, "d"]
    n = df[1, "n"]
    tex = "\\begin{tabular}{cccccccccccccccc}\n"
    tex *= "\$d=$d\$, & \$n=$n\$ &&&&&&&&&&&&\\\\\n"
    tex *= "\\hline\n"
    tex *= "\$J\$ & LS & \$B\$ & \$\\lambda\$ & RMSE & Bias & SD & Bias/SD & "
    tex *= "\$\\widehat{\\textrm{SD}}\$ & \$\\hat\\sigma^2\$ & OBias & OSD & CR & CIW \\\\\n"

    #display(df)
    for i in 1:nrow(df)
        row = df[i, :]

        if i > 1 && df[i, :J_estimator] == df[i-1, :J_estimator]
            tex *= ""
        else
            tex *= "$(df[i, :J_estimator])"
        end

        if i > 1 && df[i, :J_lifetime] == df[i-1, :J_lifetime] &&
                     df[i, :lambda_method] == df[i-1, :lambda_method]
            tex *= "&"
        else
            Jl = df[i, :J_lifetime]
            lm_fmt = lambda_method_format(df[i, :lambda_method])
            if lm_fmt != "OPT"
                hat = "\\hat"
            else
                hat = ""
            end
            tex *= "& \$$hat\\lambda_{$Jl}^{\\scriptsize{\\textrm{$lm_fmt}}}\$"
        end

        for cell in df[i, [:B_estimator, :lambda, :rmse, :bias, :sd,
                        :bias_over_sd, :sd_hat, :sigma2_hat, :bias_theory, :sd_theory,
                        :coverage, :average_width]]
            if isa(cell, Float64)
                cell = round(cell, digits=4)
            end
            tex *= "& $cell"
        end

        tex *= "\\\\\n"
    end
    tex *= "\\end{tabular}"
    write("./replication/debiasing/table_d$(d)_n$(n).tex", tex)
end

for d in unique(data[!, "d"])
    println(d)
    data_d = filter(:d => ==(d), data)
    for n in unique(data_d[!, "n"])
        println(n)
        data_d_n = filter(:n => ==(n), data_d)
        #display(data_d_n)
        make_table(data_d_n)
    end
end

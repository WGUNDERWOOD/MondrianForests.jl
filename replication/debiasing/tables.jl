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
#println(names(data))
#println(length(names(data)))
data = select!(data, sort(names(data)))
data = select!(data, [:d, :n, :J_estimator, :lambda_target, :lambda_method, :J_lifetime,
                      :B_estimator, :B_lifetime, :n_subsample,
                      :lambda_multiplier, :lambda, :rmse, :bias, :sd, :bias_over_sd,
                      :bias_theory, :sd_theory, :coverage, :average_width])
data = sort!(data, [:J_estimator, :lambda_target, order(:lambda_method, by=lambda_method_order),
                    :J_lifetime, :d, :n, :B_estimator, :B_lifetime, :n_subsample,
                    :lambda_multiplier])

#data[!, :J_lifetime] = string.(data[!, :J_lifetime])
#println(names(data))
#println(length(names(data)))
#data = sort(data, [:J_estimator, :lambda_target])
display(data)

function make_table(df)
    tex = "\\begin{tabular}{cccccccccccccc}\n"
    tex *= "\\hline\n"
    tex *= "\$J\$ & LS & \$B\$ & \$\\lambda\$ & RMSE & Bias & SD & Bias/SD & OBias & OSD & CR & CIW \\\\\n"

    display(df)
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
            # TODO format nicely
            #tex *= "& $(df[i, :J_lifetime]) $(df[i, :lambda_method])"
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
                        :bias_over_sd, :bias_theory, :sd_theory,
                        :coverage, :average_width]]
            if isa(cell, Float64)
                cell = round(cell, digits=3)
            end
            tex *= "& $cell"
        end

        tex *= "\\\\\n"
    end
    tex *= "\\end{tabular}"

    #for row in eachrow(df_Je_Jl_lm[!, [:J_estimator, :J_lifetime, :B_estimator, :lambda, :rmse, :bias, :sd,
                                       #:bias_over_sd, :bias_theory, :sd_theory,
                                       #:coverage, :average_width]])

    #for J_estimator in unique(df[!, "J_estimator"])
        #tex *= "$J_estimator "
        #df_Je = filter(:J_estimator => ==(J_estimator), df)
        #for J_lifetime in unique(df_Je[!, "J_lifetime"])
            #df_Je_Jl = filter(:J_lifetime => ==(J_lifetime), df)
            #for lambda_method in unique(df_Je_Jl[!, "lambda_method"])
                #df_Je_Jl_lm = filter(:lambda_method => ==(lambda_method), df_Je_Jl)
                ## TODO format this nicely
                #tex *= "$J_lifetime $lambda_method "
                #for row in eachrow(df_Je_Jl_lm[!, [:B_estimator, :lambda, :rmse, :bias, :sd,
                                                   #:bias_over_sd, :bias_theory, :sd_theory,
                                                   #:coverage, :average_width]])
                    #for cell in row
                        #println(cell)
                        #cell = round(cell, digits=3)
                        #tex *= "& $cell"
                    #end
                    #tex *= "\\\\\n & "
                #end
            #end
        #end
    #end

    #tex = chop(tex, tail = 3)
    #tex *= "\\hline\n"
    #tex *= "\\end{tabular}"

    n = df[1, "n"]
    d = df[1, "d"]
    write("./replication/debiasing/table_d$(d)_n$(n).tex", tex)

end

for d in unique(data[!, "d"])
    data_d = filter(:d => ==(d), data)
    for n in unique(data_d[!, "n"])
        println(n)
        data_d_n = filter(:n => ==(n), data_d)
        display(data_d_n)
        make_table(data_d_n)
    end
end

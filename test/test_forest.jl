@testset verbose = true "Forests" begin

    @testset verbose = true "Constants" begin
        n = 50
        lambda = 5.0
        n_trees = 20
        n_evals = 60
        debias_order = 0
        significance_level = 0.05
        for d in 1:3
            X = [ntuple(i -> rand(), d) for _ in 1:n]
            x_evals = [ntuple(i -> rand(), d) for _ in 1:n_evals]
            Y = [2.0 for _ in 1:n]
            forest = MondrianForest(lambda, n_trees, x_evals, debias_order,
                                    significance_level, X, Y)
            @test all(forest.mu_hat .== 2)
        end
    end


    #for d in 1:2
        #lambda = 10.0
        #n_trees = 1000
        #debias_order = 0
        #significance_level = 0.05
        #n = 2000
        #data = generate_uniform_data(d, n)
        #x_evals = [ntuple(_ -> x, d) for x in range(0.25, 0.75, length=3)]
        #forest = MondrianForest(lambda, n_trees, x_evals, debias_order,
                                #significance_level, data["X"], data["Y"])
        #show(forest)
    #end
end


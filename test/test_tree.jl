@testset verbose = true "Trees" begin
    for d in 1:3
        lambda = 2.0
        tree = MondrianTree(d, lambda)
        show(tree)
    end
end


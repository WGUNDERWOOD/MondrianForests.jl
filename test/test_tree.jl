@testset verbose = true "Trees" begin
    @testset verbose = true "Tree construction" begin
        for d in 1:3
            lambda = 3.0
            tree = MondrianTree(d, lambda)
            @test tree.creation_time == 0.0
            @test all(tree.lower .== 0.0)
            @test all(tree.upper .== 1.0)
        end
    end

    @testset verbose = true "Tree construction errors" begin
        for d in 1:3
            lambda = -1.0
            @test_throws DomainError MondrianTree(d, lambda)
        end
    end

    @testset verbose = true "is_in" begin
        for d in 1:3
            lambda = 3.0
            tree = MondrianTree(d, lambda)
            x_in = ntuple(i -> rand(), d)
            x_out = ntuple(i -> 2 + rand(), d)
            @test is_in(x_in, tree)
            @test !is_in(x_out, tree)
        end
    end

    @testset verbose = true "get_center" begin
        for d in 1:3
            lambda = 3.0
            tree = MondrianTree(d, lambda)
            @test get_center(tree) == ntuple(x -> 0.5, d)
        end
    end

    @testset verbose = true "get_volume" begin
        for d in 1:3
            lambda = 3.0
            lower = ntuple(x -> 0.3, d)
            upper = ntuple(x -> 0.7, d)
            tree = MondrianTree("", lambda, lower, upper, 0.0)
            @test isapprox(get_volume(tree), 0.4^d, rtol=1e-10)
        end
    end

    @testset verbose = true "get_common_refinement" begin
        for d in 1:3
            for B in 2:4
                lambda = 3.0
                trees = [MondrianTree(d, lambda) for b in 1:B]
                refinement = get_common_refinement(trees)
                times = [t.creation_time for tree in trees for t in get_subtrees(tree)]
                refined_times = [t.creation_time for t in get_subtrees(refinement)]
                @test Set(times) == Set(refined_times)
                @test count_leaves(refinement) >= maximum(count_leaves.(trees))
            end
        end
    end

    @testset verbose = true "are_in_same_leaf" begin
        for d in 1:3
            lambda = 10.0
            tree = MondrianTree(d, lambda)
            x_left = ntuple(i -> 0.0, d)
            x_right = ntuple(i -> 1.0, d)
            @test MondrianForests.are_in_same_leaf(x_left, x_left, tree)
            @test !MondrianForests.are_in_same_leaf(x_left, x_right, tree)
            lambda = 0.0
            tree = MondrianTree(d, lambda)
            x_left = ntuple(i -> 0.0, d)
            x_right = ntuple(i -> 1.0, d)
            @test MondrianForests.are_in_same_leaf(x_left, x_left, tree)
            @test MondrianForests.are_in_same_leaf(x_left, x_right, tree)
        end
    end

    @testset verbose = true "get_subtrees" begin
        for d in 1:3
            lambda = 3.0
            tree = MondrianTree(d, lambda)
            subtrees = get_subtrees(tree)
            @test all([t.tree_left.id[end] for t in subtrees if t.is_split] .== 'L')
            @test all([t.tree_right.id[end] for t in subtrees if t.is_split] .== 'R')
        end
    end

    @testset verbose = true "get_leaves" begin
        for d in 1:3
            lambda = 3.0
            tree = MondrianTree(d, lambda)
            leaves = get_leaves(tree)
            @test leaves == [t for t in get_subtrees(tree) if !t.is_split]
        end
    end

    @testset verbose = true "get_leaf_containing" begin
        for d in 1:3
            lambda = 3.0
            tree = MondrianTree(d, lambda)
            x = ntuple(i -> 0.5, d)
            leaf = get_leaf_containing(x, tree)
            @test is_in(x, leaf)
        end
    end

    @testset verbose = true "count_leaves" begin
        n_reps = 2000
        for d in 1:3
            for lambda in [1.0, 2.0, 3.0]
                count = 0
                for rep in 1:n_reps
                    tree = MondrianTree(d, lambda)
                    count += MondrianForests.count_leaves(tree)
                end
                count /= n_reps
                @test isapprox(count, (1 + lambda)^d, rtol=0.1)
            end
        end
    end

    @testset verbose = true "restrict" begin
        for d in 1:3
            lambda = 5.0
            tree = MondrianTree(d, lambda)
            restricted = restrict(tree, 3.0)
            times = [t.creation_time for t in get_subtrees(tree)]
            restricted_times = [t.creation_time for t in get_subtrees(restricted)]
            @test issubset(restricted_times, times)
        end
    end

    @testset verbose = true "show" begin
        for d in 1:3
            lambda = 2.0
            tree = MondrianTree(d, lambda)
            @suppress show(tree)
        end
    end
end

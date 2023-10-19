# TODO rewrite

@testset verbose = true "Trees" begin
    @testset verbose = true "Tree construction" begin
        for d in 1:3
            lambda = 3.0
            tree = MondrianTree(d, lambda)
            @test tree.creation_time == 0.0
            @test tree.cell == MondrianCell(d)
        end
    end

    @testset verbose = true "Tree construction errors" begin
        for d in 1:3
            lambda = -1.0
            @test_throws DomainError MondrianTree(d, lambda)
        end
    end

    @testset verbose = true "get_subtrees" begin
        for d in 1:3
            lambda = 3.0
            tree = MondrianTree(d, lambda)
            subtrees = get_subtrees(tree)
            println(subtrees[1])
            #@test all(collect(cell_id_left) .== 'L')
            #@test all(collect(cell_id_right) .== 'R')
        end
    end

    #=
    #
    @testset verbose = true "get_subtrees" begin
        for d in 1:3
            lambda = 2.0
            tree = MondrianTree(d, lambda)
            x_left = ntuple(i -> 0.0, d)
            x_right = ntuple(i -> 1.0, d)
            cell_id_left = MondrianForests.get_cell_id(x_left, tree)
            cell_id_right = MondrianForests.get_cell_id(x_right, tree)
            @test all(collect(cell_id_left) .== 'L')
            @test all(collect(cell_id_right) .== 'R')
        end
    end

    @testset verbose = true "are_in_same_cell" begin
        for d in 1:3
            lambda = 10.0
            tree = MondrianTree(d, lambda)
            x_left = ntuple(i -> 0.0, d)
            x_right = ntuple(i -> 1.0, d)
            @test MondrianForests.are_in_same_cell(x_left, x_left, tree)
            @test !MondrianForests.are_in_same_cell(x_left, x_right, tree)
            lambda = 0.0
            tree = MondrianTree(d, lambda)
            x_left = ntuple(i -> 0.0, d)
            x_right = ntuple(i -> 1.0, d)
            @test MondrianForests.are_in_same_cell(x_left, x_left, tree)
            @test MondrianForests.are_in_same_cell(x_left, x_right, tree)
        end
    end

    @testset verbose = true "count_cells" begin
        n_reps = 2000
        for d in 1:3
            for lambda in [1.0, 2.0, 3.0]
                count = 0
                for rep in 1:n_reps
                    tree = MondrianTree(d, lambda)
                    count += MondrianForests.count_cells(tree)
                end
                count /= n_reps
                @test isapprox(count, (1 + lambda)^d, rtol=0.1)
            end
        end
    end

    @testset verbose = true "show" begin
        for d in 1:3
            lambda = 2.0
            tree = MondrianTree(d, lambda)
            @suppress show(tree)
        end
    end
    =#
end

@testset verbose = true "Cells" begin
    @testset verbose = true "Initial cell construction" begin
        for d in 1:5
            cell = MondrianCell(d)
            @test length(cell.lower) == d
            @test length(cell.upper) == d
            @test all(cell.lower .== 0)
            @test all(cell.upper .== 1)
            lower = ntuple(i -> rand() / 2, d)
            upper = ntuple(i -> 0.5 + rand() / 2, d)
            cell = MondrianCell(lower, upper)
            @test length(cell.lower) == d
            @test length(cell.upper) == d
            @test all(0 .<= cell.lower .<= cell.upper .<= 1)
        end
    end

    @testset verbose = true "Cell construction errors" begin
        @test_throws DomainError MondrianCell(-1)
        lower = ntuple(i -> -1.0, 2)
        upper = ntuple(i -> 1.0, 2)
        @test_throws DomainError MondrianCell(lower, upper)
        lower = ntuple(i -> 1.0, 2)
        upper = ntuple(i -> -1.0, 2)
        @test_throws DomainError MondrianCell(lower, upper)
        lower = ntuple(i -> 1.0, 2)
        upper = ntuple(i -> 0.0, 2)
        @test_throws ArgumentError MondrianCell(lower, upper)
    end

    @testset verbose = true "is_in" begin
        for d in 1:5
            cell = MondrianCell(d)
            x_in = ntuple(i -> rand(), d)
            x_out = ntuple(i -> 2 + rand(), d)
            @test MondrianForests.is_in(x_in, cell)
            @test !MondrianForests.is_in(x_out, cell)
        end
    end

    @testset verbose = true "show" begin
        for d in 1:5
            lower = ntuple(i -> rand() / 2, d)
            upper = ntuple(i -> 0.5 + rand() / 2, d)
            cell = MondrianCell(lower, upper)
            @suppress show(cell)
        end
    end
end

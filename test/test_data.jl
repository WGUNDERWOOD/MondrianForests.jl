@testset verbose = true "Data" begin
    d = 2
    n = 10
    data = MondrianForests.generate_uniform_data_uniform_errors(d, n)
    data = MondrianForests.generate_uniform_data_normal_errors(d, n)
end

using MPCCModelManager, Test, SafeTestsets
@time begin
@time @safetestset "Model Procesing" begin include("test_proc_model.jl") end
end

# using MPCCModelManager
# using Test

# @testset "MPCCModelManager.jl" begin
#     # Write your tests here.
# end

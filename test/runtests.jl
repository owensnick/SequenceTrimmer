using SequenceTrimmer
using Test



@testset "SequenceTrimmer.jl" begin
    # Write your own tests here.

    @test SequenceTrimmer.count_n(dna"ACGT")  == 0
    @test SequenceTrimmer.count_n(dna"ACGTN") == 1

end

module SequenceTrimmer

using BioSequences, CodecZlib, DataStructures

export consensus, estimate_libtype, pwm_summary, trim_fastq_threads
include("trimming.jl")


end # module

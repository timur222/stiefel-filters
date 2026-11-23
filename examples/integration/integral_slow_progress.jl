include("integral.jl")

progress_interval = 20
if abspath(PROGRAM_FILE) == @__FILE__
    run_test_untegral()
end

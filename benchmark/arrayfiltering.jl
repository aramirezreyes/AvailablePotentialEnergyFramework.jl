using AvailablePotentialEnergyFramework, BenchmarkTools, Statistics


a = 100rand(512,512,80,100)
buf = similar(a)
@info "Running benchmark with nthreads = $(Threads.nthreads())"
bench =  @benchmark AvailablePotentialEnergyFramework.filter_array_2!($buf,$a,$25,$25,$1);
@info "Allocations: $(bench.allocs)"
@info "%GC time: $(bench.gctimes)"
@info "Memory: $(bench.memory)"
@info "Minimum time: $(1e-9(minimum(bench.times)))s"
@info "Median time: $(1e-9(median(bench.times)))s"
@info "Maximum time: $(1e-9(maximum(bench.times)))s"

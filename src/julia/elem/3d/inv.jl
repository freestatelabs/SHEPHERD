# Test fast matrix inversion 
using BenchmarkTools, LinearAlgebra


function inv3x3!(Ainv::AbstractArray, A::AbstractArray) 
    # Roughly 27x faster than `inv()` with no allocations

    Ainv[1,1] = (A[2,2]*A[3,3] - A[2,3]*A[3,2])/(A[1,1]*A[2,2]*A[3,3] - A[1,1]*A[2,3]*A[3,2] - A[1,2]*A[2,1]*A[3,3] + A[1,2]*A[2,3]*A[3,1] + A[1,3]*A[2,1]*A[3,2] - A[1,3]*A[2,2]*A[3,1])
    Ainv[1,2] = (-A[1,2]*A[3,3] + A[1,3]*A[3,2])/(A[1,1]*A[2,2]*A[3,3] - A[1,1]*A[2,3]*A[3,2] - A[1,2]*A[2,1]*A[3,3] + A[1,2]*A[2,3]*A[3,1] + A[1,3]*A[2,1]*A[3,2] - A[1,3]*A[2,2]*A[3,1])
    Ainv[1,3] = (A[1,2]*A[2,3] - A[1,3]*A[2,2])/(A[1,1]*A[2,2]*A[3,3] - A[1,1]*A[2,3]*A[3,2] - A[1,2]*A[2,1]*A[3,3] + A[1,2]*A[2,3]*A[3,1] + A[1,3]*A[2,1]*A[3,2] - A[1,3]*A[2,2]*A[3,1])
    Ainv[2,1] = (-A[2,1]*A[3,3] + A[2,3]*A[3,1])/(A[1,1]*A[2,2]*A[3,3] - A[1,1]*A[2,3]*A[3,2] - A[1,2]*A[2,1]*A[3,3] + A[1,2]*A[2,3]*A[3,1] + A[1,3]*A[2,1]*A[3,2] - A[1,3]*A[2,2]*A[3,1])
    Ainv[2,2] = (A[1,1]*A[3,3] - A[1,3]*A[3,1])/(A[1,1]*A[2,2]*A[3,3] - A[1,1]*A[2,3]*A[3,2] - A[1,2]*A[2,1]*A[3,3] + A[1,2]*A[2,3]*A[3,1] + A[1,3]*A[2,1]*A[3,2] - A[1,3]*A[2,2]*A[3,1])
    Ainv[2,3] = (-A[1,1]*A[2,3] + A[1,3]*A[2,1])/(A[1,1]*A[2,2]*A[3,3] - A[1,1]*A[2,3]*A[3,2] - A[1,2]*A[2,1]*A[3,3] + A[1,2]*A[2,3]*A[3,1] + A[1,3]*A[2,1]*A[3,2] - A[1,3]*A[2,2]*A[3,1])
    Ainv[3,1] = (A[2,1]*A[3,2] - A[2,2]*A[3,1])/(A[1,1]*A[2,2]*A[3,3] - A[1,1]*A[2,3]*A[3,2] - A[1,2]*A[2,1]*A[3,3] + A[1,2]*A[2,3]*A[3,1] + A[1,3]*A[2,1]*A[3,2] - A[1,3]*A[2,2]*A[3,1])
    Ainv[3,2] = (-A[1,1]*A[3,2] + A[1,2]*A[3,1])/(A[1,1]*A[2,2]*A[3,3] - A[1,1]*A[2,3]*A[3,2] - A[1,2]*A[2,1]*A[3,3] + A[1,2]*A[2,3]*A[3,1] + A[1,3]*A[2,1]*A[3,2] - A[1,3]*A[2,2]*A[3,1])
    Ainv[3,3] = (A[1,1]*A[2,2] - A[1,2]*A[2,1])/(A[1,1]*A[2,2]*A[3,3] - A[1,1]*A[2,3]*A[3,2] - A[1,2]*A[2,1]*A[3,3] + A[1,2]*A[2,3]*A[3,1] + A[1,3]*A[2,1]*A[3,2] - A[1,3]*A[2,2]*A[3,1])

end

function inv3x3_2!(Ainv::AbstractArray, A::AbstractArray) 
    # Roughly 45x faster than `inv()` with no allocations

    detA = det3x3(A)

    Ainv[1,1] = (A[2,2]*A[3,3] - A[2,3]*A[3,2])/detA
    Ainv[1,2] = (-A[1,2]*A[3,3] + A[1,3]*A[3,2])/detA
    Ainv[1,3] = (A[1,2]*A[2,3] - A[1,3]*A[2,2])/detA
    Ainv[2,1] = (-A[2,1]*A[3,3] + A[2,3]*A[3,1])/detA
    Ainv[2,2] = (A[1,1]*A[3,3] - A[1,3]*A[3,1])/detA
    Ainv[2,3] = (-A[1,1]*A[2,3] + A[1,3]*A[2,1])/detA
    Ainv[3,1] = (A[2,1]*A[3,2] - A[2,2]*A[3,1])/detA
    Ainv[3,2] = (-A[1,1]*A[3,2] + A[1,2]*A[3,1])/detA
    Ainv[3,3] = (A[1,1]*A[2,2] - A[1,2]*A[2,1])/detA

    return detA
end

function det3x3(A::AbstractArray)
    # Roughly 30x faster than `det()` with no allocations 

    return A[1,1]*A[2,2]*A[3,3] - A[1,1]*A[2,3]*A[3,2] - A[1,2]*A[2,1]*A[3,3] + A[1,2]*A[2,3]*A[3,1] + A[1,3]*A[2,1]*A[3,2] - A[1,3]*A[2,2]*A[3,1]

end

A = rand(3,3) 
Ainv = zeros(3,3) 

inv3x3!(Ainv, A) 
display(Ainv*A)

inv3x3_2!(Ainv, A)
display(Ainv*A)

println(""); println("Testing inv()"); b = @benchmark inv($A); display(b); println("")
# BenchmarkTools.Trial: 10000 samples with 205 evaluations per sample.
#  Range (min … max):  347.151 ns … 87.805 μs  ┊ GC (min … max):  0.00% … 99.30%
#  Time  (median):     461.176 ns              ┊ GC (median):     0.00%
#  Time  (mean ± σ):   722.964 ns ±  2.628 μs  ┊ GC (mean ± σ):  36.53% ± 10.62%
#   █                                                            ▁
#   █▇▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▅▆ █
#   347 ns        Histogram: log(frequency) by time      19.7 μs <
#  Memory estimate: 1.88 KiB, allocs estimate: 7.

println("Testing inv3x3!()"); b = @benchmark inv3x3!($Ainv, $A); display(b); println("")
# BenchmarkTools.Trial: 10000 samples with 998 evaluations per sample.
#  Range (min … max):  17.285 ns … 50.351 ns  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     17.368 ns              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   17.409 ns ±  0.602 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%
#        ▅     █     ▆     ▂     ▃     ▂     ▁                  ▁
#   ▅▁▁▁▁█▁▁▁▁▁█▁▁▁▁▁█▁▁▁▁▁█▁▁▁▁▁█▁▁▁▁▁█▁▁▁▁▁█▁▁▁▁▁▇▁▁▁▁▁▆▁▁▁▁▄ █
#   17.3 ns      Histogram: log(frequency) by time      17.7 ns <
#  Memory estimate: 0 bytes, allocs estimate: 0.

println("Testing inv3x3_2!()"); b = @benchmark inv3x3_2!($Ainv, $A); display(b); println("")
# BenchmarkTools.Trial: 10000 samples with 999 evaluations per sample.
#  Range (min … max):  10.761 ns … 43.043 ns  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     10.886 ns              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   11.096 ns ±  1.667 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%
#   █▆▁▄▂▁▂                                                     ▁
#   ██████████▆▄▃▃▄▃▁▁▁▃▁▁▁▁▁▃▃▁▁▁▁▄▃▁▁▃▁▃▄▄▃▁▄▄▁▁▃▃▄▁▃▁▁▁▄▃▃▄▄ █
#   10.8 ns      Histogram: log(frequency) by time      18.3 ns <
#  Memory estimate: 0 bytes, allocs estimate: 0.

println("Testing det()"); b = @benchmark det($A); display(b); println("")
# BenchmarkTools.Trial: 10000 samples with 925 evaluations per sample.
#  Range (min … max):  115.496 ns …  21.791 μs  ┊ GC (min … max):  0.00% … 99.23%
#  Time  (median):     129.639 ns               ┊ GC (median):     0.00%
#  Time  (mean ± σ):   150.358 ns ± 430.907 ns  ┊ GC (mean ± σ):  13.54% ±  5.33%
#                     ▂▃▅▅▇█▆█▇▇█▄▃▃▃▁                             
#   ▁▁▁▂▄▄▃▂▄▆▆▅▅▆▇▇████████████████████▆▆▅▅▄▄▄▄▃▃▂▂▂▂▂▂▂▂▁▁▁▁▁▁▁ ▄
#   115 ns           Histogram: frequency by time          149 ns <
#  Memory estimate: 224 bytes, allocs estimate: 4.

println("Testing det3x3()"); b = @benchmark det3x3($A); display(b); println("")
# BenchmarkTools.Trial: 10000 samples with 1000 evaluations per sample.
#  Range (min … max):  4.250 ns … 12.250 ns  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     4.375 ns              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   4.366 ns ±  0.118 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%
#            ▂         ▇        █         ▃                  ▂ ▁
#   ▄▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁█ █
#   4.25 ns      Histogram: log(frequency) by time      4.5 ns <
#  Memory estimate: 0 bytes, allocs estimate: 0.

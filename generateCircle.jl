import BenchmarkExample: Circular

for n in 3:25
    Circular.generateMsh("msh/circular_tri3_$n.msh", transfinite=n, order=1, quad=false)
# Circular.generateMsh("msh/circular_tri6_$n.msh", transfinite=n, order=2, quad=false)
end

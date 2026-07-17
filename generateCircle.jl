import BenchmarkExample: Circular

# for n in 3:25
n=0.25
    Circular.generateMsh("msh/circular_tri3_irregular_$n.msh", lc=n*2, order=1, quad=false)
    # Circular.generateMsh("msh/circular_tri3_$n.msh", transfinite=n, order=1, quad=false)
# end

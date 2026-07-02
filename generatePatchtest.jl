import BenchmarkExample: PatchTest

# PatchTest.𝐿 = 100.0
# for n in 2:25
#     PatchTest.generateMsh("msh/patchtest_tri3_$n.msh", transfinite=n+1, order=1, quad=false)
#     PatchTest.perturbMsh("msh/patchtest_tri3_$n.msh", "msh/patchtest_un_tri3_$n.msh",n ,scale=3.0, seed=n)
#     # PatchTest.generateMsh("msh/patchtest_tri6_$n.msh", transfinite=n+1, order=2, quad=false)
# # end

# for n in 2:32
#     PatchTest.generateMsh("msh/patchtest_quad4_$n.msh", transfinite=n+1, order=1, quad=true)
#     # PatchTest.generateMsh("msh/patchtest_quad8_$n.msh", transfinite=n+1, order=2, quad=true)
# end

# for n in 2:32
#     PatchTest.generateMsh("msh/patchtest_high_un_tri3_$n.msh", lc = 1/(n-1), order=1, quad=false)
# end

# for n in 3:32
#     PatchTest.generateMsh("msh/patchtest_un_quad4_$n.msh", lc = 1/(n-1), order=1, quad=true)
#     # PatchTest.generateMsh("msh/patchtest_quad8_$n.msh", transfinite=n+1, order=2, quad=true)
# end
# for n in 0.061:0.001:0.1
# for n in 2:32
# PatchTest.generateMsh("msh/patchtest_tri3_irregular_.msh", lc=0.041, order=1, quad=false)
PatchTest.generateMsh("msh/patchtest_quad4_irregular_.msh", lc=0.039, order=1, quad=true)
# end

# for n in 2:64
#     PatchTest.generateMsh("msh/patchtest_un_quad4_$n.msh", lc = 1/n, order=1, quad=true)
# end

# ---------------------------------------------------------------------------------------------------------

# for n in (0.04166)
#     PatchTest.generateMsh("msh/patchtest_high_un_tri3_$n.msh", lc = n, order=1, quad=false)
#     # PatchTest.generateMsh("msh/patchtest_tri6_$n.msh", transfinite=n+1, order=2, quad=false)
# end

# for nw in (0.1,0.2)
#     for nq in (0.3,0.4)
#     PatchTest.generateMsh("msh/patchtest_high_un_tri3_$n.msh", lc = nw*nq, order=1, quad=false)
#    # PatchTest.generateMsh("msh/patchtest_tri6_$n.msh", transfinite=n+0.001, order=2, quad=false)
#     end
# end


# n = 4
# PatchTest.generateMsh("msh/patchtest_tri3_100_$n.msh", transfinite=n+1, order=1)
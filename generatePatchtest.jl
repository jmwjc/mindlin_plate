import BenchmarkExample: PatchTest

# for n in 2:42
#     PatchTest.generateMsh("msh/patchtest_tri3_$n.msh", transfinite=n+1, order=1, quad=false)
#     PatchTest.generateMsh("msh/patchtest_tri6_$n.msh", transfinite=n+1, order=2, quad=false)
# end


for n in 2:25
    PatchTest.generateMsh("msh/patchtest_high_un_tri3_$n.msh", lc = 1/n, order=1, quad=false)
   # PatchTest.generateMsh("msh/patchtest_tri6_$n.msh", transfinite=n+0.001, order=2, quad=false)
end
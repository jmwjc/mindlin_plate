import BenchmarkExample: PatchTest

for n in 1:32
    PatchTest.generateMsh("msh/patchtest_tri3_$n.msh", transfinite=n+1, order=1, quad=false)
    PatchTest.generateMsh("msh/patchtest_tri6_$n.msh", transfinite=n+1, order=2, quad=false)
end
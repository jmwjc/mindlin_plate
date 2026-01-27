import BenchmarkExample: Circular

# 与 generatePatchtest.jl 同风格：直接调用已安装/当前环境中的 BenchmarkExample.Circular.generateMsh
# 注意：如果在默认环境运行，加载到的 BenchmarkExample 版本可能与工作区源码不同。

# tri3 / tri6
for n in 2:25
    Circular.generateMsh("msh/circular_tri3_$n.msh", transfinite=n + 1, order=1, quad=false)
    # Circular.generateMsh("msh/circular_tri6_$n.msh", transfinite=n + 1, order=2, quad=false)
end

# quad4 / quad8（可选）
# for n in 2:25
#     Circular.generateMsh("msh/circular_quad4_$n.msh", transfinite=n + 1, order=1, quad=true)
#     Circular.generateMsh("msh/circular_quad8_$n.msh", transfinite=n + 1, order=2, quad=true)
# end

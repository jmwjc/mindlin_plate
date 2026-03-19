using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ

import Gmsh: gmsh
using LinearAlgebra

# ─────────────────────────────────────────────────────────────
# patch_test_mincheck.jl
# 目的：按导师要求，用“线性/常数场”逐项验证边界条件施加是否正确。
#      不算 L2，不做自动分类，不引入复杂诊断。
#
# 几何与边界编号来自 BenchmarkExample.PatchTest.generateGeo:
#   点：(0,0)->(L,0)->(L,L)->(0,L)
#   Γ¹: (0,0)-(L,0)  => y=0  (bottom)
#   Γ²: (L,0)-(L,L)  => x=L  (right)
#   Γ³: (L,L)-(0,L)  => y=L  (top)
#   Γ⁴: (0,L)-(0,0)  => x=0  (left)
#
# 这里把它“类比圆板四分之一”的边界类型：
#   outer (类比 Γᵉ 外圆弧 clamped):  x=L  => Γ²
#   sym-x (类比 x=0 对称边，只压 βx): x=0  => Γ⁴
#   sym-y (类比 y=0 对称边，只压 βy): y=0  => Γ¹
#   其余 Γ³ 作为自然边（此脚本中不加自然边项，只用来观测是否被错误影响）
# ─────────────────────────────────────────────────────────────

# 逐项启用开关（建议一次只开一项做隔离验证）
# 仍保留 STEP 作为“默认组合方案”，但最终是否装配由下面的布尔开关决定。
const STEP = 4

# 罚项开关（优先用这些做实验隔离）
const ENFORCE_OUTER_W  = false
const ENFORCE_OUTER_PHI = false
const ENFORCE_SYMX_PHI1 = false
const ENFORCE_SYMY_PHI2 = true

# 例：实验A（只罚 Γ¹ 上 φ2）：把上面四个开关改成 false,false,false,true
# 例：实验B（只罚 Γ⁴ 上 φ1）：把上面四个开关改成 false,false,true,false

# 罚因子
const αw = 1e7
const αφ = 1e7

# 线性/常数假设场（你可以随时改成更简单，比如全 0）
# 建议：先用常数，确保“能压住”；再用线性，确保“压到正确数值”。
w_exact(x, y, z) = 2.0 + 0.3 * x - 0.7 * y

# 快速对照：交换 φ1/φ2 的 exact 定义（用于判定 d₁/d₂ 是否分量对调）
const SWAP_PHI_EXACT = false

# φ exact 模式：:constant 用于最小对照（两分量同难度），:linear 为原线性场
const PHI_MODE = :constant  # 可选 :constant / :linear

φ1_exact_linear(x, y, z) = -1.0 + 0.5 * x + 0.0 * y
φ2_exact_linear(x, y, z) = 0.2 + 0.0 * x - 0.4 * y

φ1_exact_const(x, y, z) = 1.0
φ2_exact_const(x, y, z) = -2.0

φ1_exact_base(x, y, z) = (PHI_MODE === :constant) ? φ1_exact_const(x, y, z) : φ1_exact_linear(x, y, z)
φ2_exact_base(x, y, z) = (PHI_MODE === :constant) ? φ2_exact_const(x, y, z) : φ2_exact_linear(x, y, z)

φ1_exact(x, y, z) = SWAP_PHI_EXACT ? φ2_exact_base(x, y, z) : φ1_exact_base(x, y, z)
φ2_exact(x, y, z) = SWAP_PHI_EXACT ? φ1_exact_base(x, y, z) : φ2_exact_base(x, y, z)

# 域内载荷：排查边界施加时，建议先置零，避免域内项拉扯边界
q_zero(x, y, z) = 0.0
m1_zero(x, y, z) = 0.0
m2_zero(x, y, z) = 0.0

E = 1.0
ν = 0.3
h = 1e-1

# ─────────────────────────────────────────────────────────────
# 工具：边界最大残差（只看唯一节点）
# ─────────────────────────────────────────────
function max_res_on_boundary(elements_Γ, nodes, getter, exact_fun)
    maxv = 0.0
    seen = Dict{Int,Bool}()
    for elm in elements_Γ
        for nd in elm.𝓒
            if get(seen, nd.𝐼, false)
                continue
            end
            seen[nd.𝐼] = true
            maxv = max(maxv, abs(getter(nd) - exact_fun(nd.x, nd.y, nd.z)))
        end
    end
    return maxv
end

# ─────────────────────────────────────────────────────────────
# 主程序
# ─────────────────────────────────────────────────────────────

gmsh.initialize()
gmsh.open("msh/patchtest.msh")
entities = getPhysicalGroups()
nodes = get𝑿ᵢ()

n = length(nodes)

kww = zeros(n, n)
kqq = zeros(2 * n, 2 * n)
kqw = zeros(2 * n, n)
fw = zeros(n)
fq = zeros(2 * n)

# 域内组装
elementsΩ = getElements(nodes, entities["Ω"])
prescribe!(elementsΩ, :E => E, :ν => ν, :h => h, :q => q_zero, :m₁ => m1_zero, :m₂ => m2_zero)
set∇𝝭!(elementsΩ)

(∫wwdΩ => elementsΩ)(kww)
(∫φwdΩ => elementsΩ)(kqw)
((∫φφdΩ => elementsΩ))(kqq)
((∫κκdΩ => elementsΩ))(kqq)

(∫wqdΩ => elementsΩ)(fw)
(∫φmdΩ => elementsΩ)(fq)

# 边界 elements
Γ1 = getElements(nodes, entities["Γ¹"])
Γ2 = getElements(nodes, entities["Γ²"])
Γ3 = getElements(nodes, entities["Γ³"])
Γ4 = getElements(nodes, entities["Γ⁴"])
set𝝭!(Γ1);
set𝝭!(Γ2);
set𝝭!(Γ3);
set𝝭!(Γ4);

# outer/sym 类比圆板 clamped
Γ_outer = Γ2   # x=L
Γ_symx = Γ4   # x=0
Γ_symy = Γ1   # y=0

println("[map] outer=Γ² (x=L), symx=Γ⁴ (x=0), symy=Γ¹ (y=0), natural=Γ³ (y=L)")
println("[param] STEP=", STEP, ", αw=", αw, ", αφ=", αφ)

# 边界最小诊断：元素数和权重总和（确认该边界确实有积分贡献）
function boundary_stats(elements_Γ)
    ne = length(elements_Γ)
    wsum = 0.0
    for elm in elements_Γ
        for ξ in elm.𝓖
            wsum += ξ.𝑤
        end
    end
    return ne, wsum
end

let
    ne2, w2 = boundary_stats(Γ_outer)
    ne4, w4 = boundary_stats(Γ_symx)
    ne1, w1 = boundary_stats(Γ_symy)
    println("[diag-bnd] Γ² outer: ne=", ne2, ", Σw=", w2)
    println("[diag-bnd] Γ⁴ symx : ne=", ne4, ", Σw=", w4)
    println("[diag-bnd] Γ¹ symy : ne=", ne1, ", Σw=", w1)
end

# outer 上 w
if ENFORCE_OUTER_W
    kww0 = copy(kww); fw0 = copy(fw)
    prescribe!(Γ_outer, :α => αw * E, :g => w_exact)
    (∫αwwdΓ => Γ_outer)(kww, fw)
    println("[diag-asm] Γ² ∫αwwdΓ: Δ‖kww‖=", norm(kww - kww0), ", Δ‖fw‖=", norm(fw - fw0))
end

# outer 上 φ1/φ2
if ENFORCE_OUTER_PHI
    kqq0 = copy(kqq); fq0 = copy(fq)
    prescribe!(Γ_outer, :α => αφ * E, :g₁ => φ1_exact, :g₂ => φ2_exact, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
    (∫αφφdΓ => Γ_outer)(kqq, fq)
    println("[diag-asm] Γ² ∫αφφdΓ: Δ‖kqq‖=", norm(kqq - kqq0), ", Δ‖fq‖=", norm(fq - fq0))
end

# sym-x 上只压 φ1
if ENFORCE_SYMX_PHI1
    kqq0 = copy(kqq); fq0 = copy(fq)
    prescribe!(Γ_symx, :α => αφ * E, :g₁ => φ1_exact, :g₂ => φ2_exact, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 0.0)
    (∫αφφdΓ => Γ_symx)(kqq, fq)
    println("[diag-asm] Γ⁴ ∫αφφdΓ: Δ‖kqq‖=", norm(kqq - kqq0), ", Δ‖fq‖=", norm(fq - fq0))
end

# sym-y 上只压 φ2
if ENFORCE_SYMY_PHI2
    kqq0 = copy(kqq); fq0 = copy(fq)
    prescribe!(Γ_symy, :α => αφ * E, :g₁ => φ1_exact, :g₂ => φ2_exact, :n₁₁ => 0.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
    (∫αφφdΓ => Γ_symy)(kqq, fq)
    println("[diag-asm] Γ¹ ∫αφφdΓ: Δ‖kqq‖=", norm(kqq - kqq0), ", Δ‖fq‖=", norm(fq - fq0))
end

# 求解并回填
K = [kqq kqw; kqw' kww]
F = [fq; fw]

# 诊断：K 是否近奇异/是否明显非对称/是否出现 NaN/Inf
sym_err = norm(K - K') / max(norm(K), eps())
has_nan = any(isnan, K) || any(isnan, F)
has_inf = any(isinf, K) || any(isinf, F)
println("\n[diag] K symmetry rel. error = ", sym_err)
println("[diag] hasNaN(K/F)          = ", has_nan)
println("[diag] hasInf(K/F)          = ", has_inf)

# rcond 有时会因数值尺度/分解过程返回 NaN；这里改用奇异值估计 cond
let
    try
        σ = svdvals(K)
        σmax = maximum(σ)
        σmin = minimum(σ)
        cond_est = σmax / σmin
        println("[diag] σmin(K)              = ", σmin)
        println("[diag] σmax(K)              = ", σmax)
        println("[diag] cond_est(K)=σmax/σmin = ", cond_est)
    catch e
        println("[diag] svdvals(K) failed: ", typeof(e), " / ", e)
        println("[diag] σmin(K)              = NaN")
        println("[diag] σmax(K)              = NaN")
        println("[diag] cond_est(K)=σmax/σmin = NaN")
    end
end
println("[diag] max|F|                = ", maximum(abs.(F)))

sol = K \ F
println("[diag] max|sol|             = ", maximum(abs.(sol)))

push!(nodes, :d => sol[2 * n + 1:end], :d₁ => sol[1:2:2 * n], :d₂ => sol[2:2:2 * n])

# 输出：只看被施加的边界残差
# 注意：若 STEP<某项，该边界残差输出就不代表该边界被约束。
println("\n[check] boundary max residuals:")
if ENFORCE_OUTER_W
    println("  outer Γ²: max|w-w_exact|   = ", max_res_on_boundary(Γ_outer, nodes, nd -> nd.d, w_exact))
end
if ENFORCE_OUTER_PHI
    println("  outer Γ²: max|φ1-φ1_exact| = ", max_res_on_boundary(Γ_outer, nodes, nd -> nd.d₁, φ1_exact))
    println("  outer Γ²: max|φ2-φ2_exact| = ", max_res_on_boundary(Γ_outer, nodes, nd -> nd.d₂, φ2_exact))
end
if ENFORCE_SYMX_PHI1
    println("  symx  Γ⁴: max|φ1-φ1_exact| = ", max_res_on_boundary(Γ_symx, nodes, nd -> nd.d₁, φ1_exact))
end
if ENFORCE_SYMY_PHI2
    println("  symy  Γ¹: max|φ2-φ2_exact| = ", max_res_on_boundary(Γ_symy, nodes, nd -> nd.d₂, φ2_exact))
end

# observe: 在“只压单分量”时，观察另一个分量是否被误压（应当不必然接近 exact）
if ENFORCE_SYMX_PHI1 && !ENFORCE_SYMY_PHI2
    println("  [observe] symx Γ⁴: max|φ2-φ2_exact| = ", max_res_on_boundary(Γ_symx, nodes, nd -> nd.d₂, φ2_exact))
end
if ENFORCE_SYMY_PHI2 && !ENFORCE_SYMX_PHI1
    println("  [observe] symy Γ¹: max|φ1-φ1_exact| = ", max_res_on_boundary(Γ_symy, nodes, nd -> nd.d₁, φ1_exact))
end

# 额外：自然边（Γ³）不应被你“误加”强制项（这里仅观测值，不判定对错）
println("\n[observe] natural Γ³ values (not enforced in this script):")
println("  Γ³: max|w|  = ", maximum(abs.([nd.d for elm in Γ3 for nd in elm.𝓒])))
println("  Γ³: max|φ1| = ", maximum(abs.([nd.d₁ for elm in Γ3 for nd in elm.𝓒])))
println("  Γ³: max|φ2| = ", maximum(abs.([nd.d₂ for elm in Γ3 for nd in elm.𝓒])))

gmsh.finalize()

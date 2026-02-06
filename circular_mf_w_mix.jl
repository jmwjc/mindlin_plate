using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QQdΩ, ∫Q∇wdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ

using LinearAlgebra: dot
using TimerOutputs, XLSX
import Gmsh: gmsh

# ------------------------------------------------------------
# 圆板算例（Katili 1993）
# 工况：1/4 圆板，R=5；中厚板 h=1（R/h=5）；外圆弧固接；均布载荷 fz=1
# msh: msh/circular_tri3_16.msh（必须包含 Ω, Γ, Γᵇ, Γᵉ, Γˡ, 𝐴）
# ------------------------------------------------------------

E = 10.92
ν = 0.3
h = 1.0
R = 5.0
kappa = 5 / 6
Dᵇ = E * h^3 / 12 / (1 - ν^2)
Dˢ = kappa * E * h / (2 * (1 + ν))

# 载荷（均布）
F = 1.0
q(x, y, z) = F

# 论文给出的中心位移参考值（Table VIa: clamped, R/h=5）
w_center_exact = 11.551
M_center_exact = 2.03
U_total_exact = 81.45

# Dirichlet 目标函数
w_bc(x, y, z) = 0.0
φ1_bc(x, y, z) = 0.0
φ2_bc(x, y, z) = 0.0

const to = TimerOutput()

gmsh.initialize()

integrationOrder = 4

type_w = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_φ = :tri3
type_q = :(PiecewisePolynomial{:Quadratic2D})
ndiv = 9
# ndiv_w = 9
XLSX.openxlsx("xls/circular_mf_w_mix.xlsx", mode="w") do xf
    for ndiv_w = 9:16
        row = ndiv_w
        # ──────────────────────────────────────────────────────────
        @timeit to "open msh file" gmsh.open("msh/circular_tri3_$ndiv_w.msh")
        @timeit to "get nodes" nodes_w = get𝑿ᵢ()
        xʷ = nodes_w.x
        yʷ = nodes_w.y
        zʷ = nodes_w.z
        nʷ = length(nodes_w)
        # RK 邻域：圆板扇形在圆弧/角点附近更易邻域退化，导致 moment matrix 数值非 SPD。
        # 不修改 ApproxOperator 的前提下，只能在脚本侧提高邻域稳定性。
        sp = RegularGrid(xʷ, yʷ, zʷ, n=1, γ=2)
        s = 1 / ndiv_w
        s₁ = 2 * s * ones(nʷ)
        s₂ = 2 * s * ones(nʷ)
        s₃ = 2 * s * ones(nʷ)
        push!(nodes_w, :s₁ => s₁, :s₂ => s₂, :s₃ => s₃)


        @timeit to "open msh file" gmsh.open("msh/circular_tri3_$ndiv.msh")
        @timeit to "get nodes" nodes_φ = get𝑿ᵢ()
        @timeit to "get entities" entities = getPhysicalGroups()

        # 关键物理组检查
        for key in ("Ω", "Γ", "Γᵇ", "Γᵉ", "Γˡ", "𝐴")
            haskey(entities, key) || error("Mesh physical group '$key' not found in msh/circular_tri3_$ndiv_w.msh")
        end

        @timeit to "calculate main elements" begin
            @timeit to "get elements" elements_φ = getElements(nodes_φ, entities["Ω"], integrationOrder)
            @timeit to "get elements" elements_w = getElements(nodes_w, entities["Ω"], eval(type_w), integrationOrder, sp)
            @timeit to "get elements" elements_q = getPiecewiseElements(entities["Ω"], eval(type_q), integrationOrder)
        end
        nₑ = length(elements_φ)
        nᵠ = length(nodes_φ)
        nᵛ = nₑ * ApproxOperator.get𝑛𝑝(elements_q[1])
        kʷʷ = zeros(nʷ, nʷ)
        kᵛʷ = zeros(2 * nᵛ, nʷ)
        kᵠʷ = zeros(2 * nᵠ, nʷ)
        fʷ = zeros(nʷ)
        kᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)
        kᵛᵛ = zeros(2 * nᵛ, 2 * nᵛ)
        kᵛᵠ = zeros(2 * nᵛ, 2 * nᵠ)
        fᵠ = zeros(2 * nᵠ)
        fᵛ = zeros(2 * nᵛ)

        # ------------------------------------------------------------
        # 组装：域积分 + 混合边界耦合（保持模板结构）
        # ------------------------------------------------------------
        @timeit to "calculate domain" begin
            @timeit to "get elements" elements_w_Γ = getElements(nodes_w, entities["Γ"], eval(type_w), integrationOrder, sp, normal=true)
            @timeit to "get elements" elements_q_Γ = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], eval(type_q), integrationOrder)

            prescribe!(elements_φ, :E => E, :ν => ν, :h => h)
            prescribe!(elements_q, :E => E, :ν => ν, :h => h)
            prescribe!(elements_w, :E => E, :ν => ν, :h => h, :q => q)
            # 圆板算例无分布弯矩载荷，显式设为 0，避免 ∫φmdΩ 读取不到 m₁/m₂
            prescribe!(elements_φ, :m₁ => (x, y, z) -> 0.0, :m₂ => (x, y, z) -> 0.0)

            @timeit to "calculate shape functions" set∇𝝭!(elements_φ)
            @timeit to "calculate shape functions" set∇𝝭!(elements_q)
            @timeit to "calculate shape functions" set∇𝝭!(elements_w)
            @timeit to "calculate shape functions" set𝝭!(elements_w_Γ)
            @timeit to "calculate shape functions" set𝝭!(elements_q_Γ)

            𝑎ᵠᵠ = ∫κκdΩ => elements_φ
            𝑎ᵛᵠ = ∫QφdΩ => (elements_q, elements_φ)
            𝑎ᵛᵛ = ∫QQdΩ => elements_q
            𝑎ᵛʷ = [
                ∫∇QwdΩ => (elements_q, elements_w),
                ∫QwdΓ => (elements_q_Γ, elements_w_Γ),
            ]
            𝑓ᵠ_op = ∫φmdΩ => elements_φ
            𝑓ʷ_op = ∫wqdΩ => elements_w

            @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)
            @timeit to "assemble" 𝑎ᵛᵛ(kᵛᵛ)
            @timeit to "assemble" 𝑎ᵛᵠ(kᵛᵠ)
            @timeit to "assemble" 𝑎ᵛʷ(kᵛʷ)
            @timeit to "assemble" 𝑓ᵠ_op(fᵠ)
            @timeit to "assemble" 𝑓ʷ_op(fʷ)
        end

        # ------------------------------------------------------------
        # 边界条件（论文）：
        # - 外圆弧 Γᵉ：固接 w=0, φ₁=0, φ₂=0
        # - 对称边界：on CB: βy=0；on CA: βx=0
        #   结合本 msh：Γˡ(x=0) -> 约束 φ₂=0；Γᵇ(y=0) -> 约束 φ₁=0
        # ------------------------------------------------------------
        α = 1e8 * E

        @timeit to "penalty on Γᵉ (w)" begin
            elements_w_Γe = getElements(nodes_w, entities["Γᵉ"], eval(type_w), integrationOrder, sp, normal=true)
            prescribe!(elements_w_Γe, :α => α, :g => w_bc)
            set𝝭!(elements_w_Γe)
            𝑎w = ∫αwwdΓ => elements_w_Γe
            𝑎w(kʷʷ, fʷ)
        end

        @timeit to "penalty on Γᵉ (φ)" begin
            elements_φ_Γe = getElements(nodes_φ, entities["Γᵉ"], integrationOrder, normal=true)
            prescribe!(elements_φ_Γe, :α => α, :g => w_bc, :g₁ => φ1_bc, :g₂ => φ2_bc, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
            set𝝭!(elements_φ_Γe)
            𝑎φ = ∫αφφdΓ => elements_φ_Γe
            𝑎φ(kᵠᵠ, fᵠ)
        end

        @timeit to "symmetry on Γᵇ (φ₁=0)" begin
            elements_φ_Γb = getElements(nodes_φ, entities["Γᵇ"], integrationOrder, normal=true)
            # Γᵇ: y=0，对称边，约束 φ₁
            prescribe!(elements_φ_Γb, :α => α, :g => w_bc, :g₁ => φ1_bc, :g₂ => w_bc, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 0.0)
            set𝝭!(elements_φ_Γb)
            𝑎φb = ∫αφφdΓ => elements_φ_Γb
            𝑎φb(kᵠᵠ, fᵠ)
        end

        @timeit to "symmetry on Γˡ (φ₂=0)" begin
            elements_φ_Γl = getElements(nodes_φ, entities["Γˡ"], integrationOrder, normal=true)
            # Γˡ: x=0，对称边，约束 φ₂
            prescribe!(elements_φ_Γl, :α => α, :g => w_bc, :g₁ => w_bc, :g₂ => φ2_bc, :n₁₁ => 0.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
            set𝝭!(elements_φ_Γl)
            𝑎φl = ∫αφφdΓ => elements_φ_Γl
            𝑎φl(kᵠᵠ, fᵠ)
        end

        # ------------------------------------------------------------
        # 求解
        # ------------------------------------------------------------
        @timeit to "solve" d = [kᵠᵠ kᵠʷ kᵛᵠ'; kᵠʷ' kʷʷ kᵛʷ'; kᵛᵠ kᵛʷ kᵛᵛ] \ [fᵠ; fʷ; fᵛ]

        nodes_q = 𝑿ᵢ[]
        for elm in elements_q
            for node in elm.𝓒
                push!(nodes_q, node)
            end
        end

        push!(nodes_φ, :d₁ => d[1:2:2*nᵠ], :d₂ => d[2:2:2*nᵠ])
        push!(nodes_w, :d => d[2*nᵠ+1:2*nᵠ+nʷ])
        push!(nodes_q, :q₁ => d[2*nᵠ+nʷ+1:2:end], :q₂ => d[2*nᵠ+nʷ+2:2:end])

        # ------------------------------------------------------------
        # 收敛指标（Katili 1993, clamped）
        # - 中心位移: w_center
        # - 中心弯矩: M_center（取中心节点对应的板弯矩）
        # - 总能量: U_total = 0.5 * f^T d（线性系统的外力做功）
        # ------------------------------------------------------------
        I_center = 1
        w_center = nodes_w.d[I_center]
        err_w = abs(w_center - w_center_exact) / abs(w_center_exact)

        # 由中心节点近似中心弯矩（如果你有更精确的中心节点定位，可替换 I_center）
        # 这里使用各向同性板的相关量从元素场恢复通常需要额外后处理；
        # 优先尝试从 nodes_φ 上的转角二阶导数恢复弯矩并不直接可用。
        # 因此此处先留默认值 NaN，避免误导；你如果已在 ApproxOperator 中有现成算子可直接替换。
        M_center = NaN
        err_M = (isfinite(M_center) && M_center_exact != 0) ? abs(M_center - M_center_exact) / abs(M_center_exact) : NaN

        # 总能量（外力做功）
        f_all = [fᵠ; fʷ; fᵛ]
        U_total = 0.5 * dot(f_all, d)
        err_U = abs(U_total - U_total_exact) / abs(U_total_exact)

        println(to)
        println("center w: ", w_center, ", exact: ", w_center_exact, ", rel err: ", err_w)
        println("total energy: ", U_total, ", exact: ", U_total_exact, ", rel err: ", err_U)
        println("center moment M (placeholder): ", M_center, ", exact: ", M_center_exact)

        # ──────────────────────────────────────────────────────────
        sheet = xf[1]
        XLSX.rename!(sheet, "circular")
        sheet["A1"] = "type w"
        sheet["B1"] = "nʷ"
        sheet["C1"] = "type φ"
        sheet["D1"] = "nᵠ"
        sheet["E1"] = "type Q"
        sheet["F1"] = "nᵛ"
        sheet["G1"] = "rel_err_w_center"
        sheet["H1"] = "rel_err_M_center"
        sheet["I1"] = "rel_err_U_total"

        # sheet["J1"] = "w_center"
        # sheet["K1"] = "w_exact"
        # sheet["L1"] = "M_center"
        # sheet["M1"] = "M_exact"
        # sheet["N1"] = "U_total"
        # sheet["O1"] = "U_exact"

        sheet["A$row"] = "$type_w"
        sheet["B$row"] = nʷ
        sheet["C$row"] = "$type_φ"
        sheet["D$row"] = nᵠ
        sheet["E$row"] = "$type_q"
        sheet["F$row"] = nᵛ
        sheet["G$row"] = log10(err_w)
        sheet["H$row"] = (isfinite(err_M) ? log10(err_M) : "")
        sheet["I$row"] = log10(err_U)

        # sheet["J$row"] = w_center
        # sheet["K$row"] = w_center_exact
        # sheet["L$row"] = M_center
        # sheet["M$row"] = M_center_exact
        # sheet["N$row"] = U_total
        # sheet["O$row"] = U_total_exact
    end
end

gmsh.finalize()

using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, L₂, L₂φ, L₂Q

using TimerOutputs, XLSX
import Gmsh: gmsh

# ─────────────────────────────────────────────────────────────
# 圆形板（四分之一圆）— 固支(clamped) — 中厚板（R/h=5）
# 参考：/home/jason/vscode/Group (自用版)/Documents/圆板精確解.md
# 物理组来自：mindlin_plate/msh/circular.geo
#   Ω: 面域
#   Γ: 全边界（合并）
#   Γᵇ: bottom 径向边（x 轴）
#   Γˡ: left   径向边（y 轴）
#   Γᵉ: 外圆弧边界（r=R）
#           𝐴: 圆心点
# ─────────────────────────────────────────────────────────────



# 几何/材料/载荷
R = 5.0
E = 10.92
ν = 0.3
k = 5 / 6
h = 1.0            # 中厚板：R/h = 5
f_z = 1.0

Dᵇ = E * h^3 / 12 / (1 - ν^2)
Dˢ = k * E * h / (2 * (1 + ν))

# 无量纲半径
ζ(x, y) = (x^2 + y^2)^0.5 / R

# ─────────────────────────────────────────────────────────────
# clamped 精确解（文档给的是轴对称 w(r)；这里把 r->sqrt(x^2+y^2)）
# 注意：该精确解对应“整圆”均布载荷；四分之一圆 + 对称边界应与之相容。
# w_exact 用于误差 L₂；φ_exact 用于 L₂φ；Q_exact 用于 L₂Q。
#
# 对于轴对称：φ 为转角向量（φx, φy）= (dw/dr) * (x/r, y/r)
# Q 为横向剪力向量 = (T_r * x/r, T_r * y/r)
# ─────────────────────────────────────────────────────────────

# w(r)
function w_exact(x, y, z)
    r = (x^2 + y^2)^0.5
    ζv = r / R
    return (f_z * R^4) / (64 * Dᵇ) * (1 - ζv^2) * ((1 - ζv^2) + 8 * (h / R)^2 / (3 * k * (1 - ν)))
end

# dw/dr
function dwdr_exact(r)
    # w = C*(1-ζ^2)*((1-ζ^2)+A)
    # with ζ=r/R, C=fz*R^4/(64Db), A=8*(h/R)^2/(3k(1-ν))
    C = (f_z * R^4) / (64 * Dᵇ)
    A = 8 * (h / R)^2 / (3 * k * (1 - ν))
    ζv = r / R
    # Let u = (1-ζ^2); w = C*u*(u + A) = C*(u^2 + A*u)
    # dw/dζ = C*(2*u*du/dζ + A*du/dζ) = C*(2*u + A)*du/dζ
    # du/dζ = -2ζ
    # dw/dr = (dw/dζ)*(dζ/dr) = C*(2*u + A)*(-2ζ)*(1/R)
    u = 1 - ζv^2
    return C * (2 * u + A) * (-2 * ζv) / R
end

function φ_exact_components(x, y)
    r = (x^2 + y^2)^0.5
    if r == 0.0
        return 0.0, 0.0
    end
    t = dwdr_exact(r)
    return t * (x / r), t * (y / r)
end

φ₁_exact(x, y, z) = (φ_exact_components(x, y))[1]
φ₂_exact(x, y, z) = (φ_exact_components(x, y))[2]

# T_r = -f_z*r/2
function Q_exact_components(x, y)
    r = (x^2 + y^2)^0.5
    if r == 0.0
        return 0.0, 0.0
    end
    Tr = -f_z * r / 2
    return Tr * (x / r), Tr * (y / r)
end

Q₁_exact(x, y, z) = (Q_exact_components(x, y))[1]
Q₂_exact(x, y, z) = (Q_exact_components(x, y))[2]

# 本圆板算例：采用体载形式 q=f_z（与模板的 patch test 从解析反推不同）
q_load(x, y, z) = f_z

# 固支边界值：圆周 r=R 上 w=0, φ=0
w_bc(x, y, z) = 0.0
φ₁_bc(x, y, z) = 0.0
φ₂_bc(x, y, z) = 0.0

# ─────────────────────────────────────────────────────────────
# 离散设置（保持模板风格）
# ─────────────────────────────────────────────────────────────
const to = TimerOutput()
gmsh.initialize()

integrationOrder = 4
integrationOrder_err = 10

type_w = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_φ = :tri3
type_q = :(PiecewisePolynomial{:Quadratic2D})
# φ/Q 使用固定网格；w 用可变网格（保持与你的模板一致）
ndiv_φ = 15
# ndiv_w = 15
XLSX.openxlsx("xls/circular_clamped_mid.xlsx", mode="w") do xf
     for ndiv_w in (3,6,9,12,15)
    row = ndiv_w

    # 1) w 网格（RK）
    @timeit to "open w msh" gmsh.open("msh/circular_tri3_$ndiv_w.msh")
    @timeit to "get nodes_w" nodes_w = get𝑿ᵢ()
    xʷ = nodes_w.x
    yʷ = nodes_w.y
    zʷ = nodes_w.z
    nʷ = length(nodes_w)

    sp = RegularGrid(xʷ, yʷ, zʷ, n=3, γ=5)
    s = (R) / (ndiv_w) # 以半径尺度做一个相对一致的支撑域
    s₁ = 1.5 * s * ones(nʷ)
    s₂ = 1.5 * s * ones(nʷ)
    s₃ = 1.5 * s * ones(nʷ)
    push!(nodes_w, :s₁ => s₁, :s₂ => s₂, :s₃ => s₃)

    # 2) φ/Q 网格（固定）
    @timeit to "open φ msh" gmsh.open("msh/circular_tri3_$ndiv_φ.msh")
    @timeit to "get nodes_φ" nodes_φ = get𝑿ᵢ()
    @timeit to "get entities" entities = getPhysicalGroups()

    # 物理组 Γ 在 msh 里可能与 Γᵇ/Γᵉ/Γˡ 重叠并导致边界元素重复；这里显式用三段边界并集代替 Γ。
    Γall = let
        dim_b, tags_b = entities["Γᵇ"]
        dim_e, tags_e = entities["Γᵉ"]
        dim_l, tags_l = entities["Γˡ"]
        @assert dim_b == 1 && dim_e == 1 && dim_l == 1
        1 => unique(vcat(tags_b, tags_e, tags_l))
    end

    # 3) 主域 elements
    @timeit to "calculate main elements" begin
        elements_φ = getElements(nodes_φ, entities["Ω"], integrationOrder)
        elements_w = getElements(nodes_w, entities["Ω"], eval(type_w), integrationOrder, sp)
        elements_q = getPiecewiseElements(entities["Ω"], eval(type_q), integrationOrder)
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

    @timeit to "calculate ∫κκdΩ" begin
        @timeit to "get elements" elements_w_Γ = getElements(nodes_w, Γall, eval(type_w), integrationOrder, sp, normal=true)
        @timeit to "get elements" elements_q_Γ = getPiecewiseBoundaryElements(Γall, entities["Ω"], eval(type_q), integrationOrder)

        prescribe!(elements_φ, :E => E, :ν => ν, :h => h)
        prescribe!(elements_q, :E => E, :ν => ν, :h => h)
        prescribe!(elements_w, :E => E, :ν => ν, :h => h, :q => q_load)

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
        𝑓ʷ_ = ∫wqdΩ => elements_w

        @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)
        @timeit to "assemble" 𝑎ᵛᵛ(kᵛᵛ)
        @timeit to "assemble" 𝑎ᵛᵠ(kᵛᵠ)
        @timeit to "assemble" 𝑎ᵛʷ(kᵛʷ)
        @timeit to "assemble" 𝑓ʷ_(fʷ)
    end

    @timeit to "calculate ∫αwwdΓ" begin
        @timeit to "get elements" elements_w_bc = getElements(nodes_w, Γall, eval(type_w), integrationOrder, sp, normal=true)
        prescribe!(elements_w_bc, :α => 1e8 * E, :g => w_bc)
        @timeit to "calculate shape functions" set𝝭!(elements_w_bc)
        𝑎ʷ = ∫αwwdΓ => elements_w_bc
        @timeit to "assemble" 𝑎ʷ(kʷʷ, fʷ)
    end

    @timeit to "calculate ∫αφφdΓ" begin
        @timeit to "get elements" elements_φ_bc = getElements(nodes_φ, Γall, integrationOrder, normal=true)
        prescribe!(elements_φ_bc, :α => 1e8 * E, :g₁ => φ₁_bc, :g₂ => φ₂_bc, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
        @timeit to "calculate shape functions" set𝝭!(elements_φ_bc)
        𝑎ᵠ = ∫αφφdΓ => elements_φ_bc
        @timeit to "assemble" 𝑎ᵠ(kᵠᵠ, fᵠ)
    end

    @timeit to "calculate ∫QwdΓ" begin
        @timeit to "get elements" elements_q_Γ = getPiecewiseBoundaryElements(Γall, entities["Ω"], eval(type_q), integrationOrder)
        # elements_q_Γ 已由 Γall 生成，无需再用 Γall 自筛选一次（该函数依赖元素拼接顺序，易出现越界）
        @timeit to "get elements" elements_w_bc = getElements(nodes_w, Γall, eval(type_w), integrationOrder, sp, normal=true)

        prescribe!(elements_w_bc, :g => w_bc)
        @timeit to "calculate shape functions" set𝝭!(elements_q_Γ)
        @timeit to "calculate shape functions" set𝝭!(elements_w_bc)

        𝑎ᵛ = ∫QwdΓ => (elements_q_Γ, elements_w_bc)
        @timeit to "assemble" 𝑎ᵛ(kᵛʷ, fᵛ)
    end

    # 6) 求解
    @timeit to "solve" d = [kᵠᵠ kᵠʷ kᵛᵠ';
        kᵠʷ' kʷʷ kᵛʷ';
        kᵛᵠ kᵛʷ kᵛᵛ] \ [fᵠ; fʷ; fᵛ]

    # 7) 回填（模板同款）
    nodes_q = 𝑿ᵢ[]
    for elm in elements_q
        for node in elm.𝓒
            push!(nodes_q, node)
        end
    end

    push!(nodes_φ, :d₁ => d[1:2:2*nᵠ], :d₂ => d[2:2:2*nᵠ])
    push!(nodes_w, :d => d[2*nᵠ+1:2*nᵠ+nʷ])
    push!(nodes_q, :q₁ => d[2*nᵠ+nʷ+1:2:end], :q₂ => d[2*nᵠ+nʷ+2:2:end])

    # 8) 误差（用精确解）
    @timeit to "calculate error" begin
        @timeit to "get elements" elements_φ_err = getElements(nodes_φ, entities["Ω"], integrationOrder_err)
        @timeit to "get elements" elements_w_err = getElements(nodes_w, entities["Ω"], eval(type_w), integrationOrder_err, sp)
        @timeit to "get elements" elements_q_err = getPiecewiseElements(entities["Ω"], eval(type_q), integrationOrder_err)

        # teacher's workaround: L₂Q 需要节点上存在 q₁/q₂；误差算例不关注该项时，补齐为 0 以避免 KeyError
        # 注意：对 Vector{Node} 调用 push! 只会写入第一个节点；这里按节点逐个补齐。
        uniq_nodes_q_err = Dict{Int,𝑿ᵢ}()
        for elm in elements_q_err
            for node in elm.𝓒
                uniq_nodes_q_err[node.𝐼] = node
            end
        end
        n_q_err = isempty(uniq_nodes_q_err) ? 0 : maximum(keys(uniq_nodes_q_err))
        q₁_err = zeros(n_q_err)
        q₂_err = zeros(n_q_err)
        for node in values(uniq_nodes_q_err)
            push!(node, :q₁ => q₁_err, :q₂ => q₂_err)
        end

        prescribe!(elements_φ_err, :E => E, :ν => ν, :h => h, :φ₁ => φ₁_exact, :φ₂ => φ₂_exact)
        prescribe!(elements_w_err, :E => E, :ν => ν, :h => h, :u => w_exact)
        prescribe!(elements_q_err, :E => E, :ν => ν, :h => h, :Q₁ => Q₁_exact, :Q₂ => Q₂_exact)

        @timeit to "calculate shape functions" set𝝭!(elements_φ_err)
        @timeit to "calculate shape functions" set𝝭!(elements_w_err)
        @timeit to "calculate shape functions" set𝝭!(elements_q_err)

        L₂_w = L₂(elements_w_err)
        L₂_φ = L₂φ(elements_φ_err)
        L₂_Q = L₂Q(elements_q_err)
    end

    println(to)

    println("L₂ error of w: ", L₂_w)
    println("L₂ error of φ: ", L₂_φ)
    println("L₂ error of Q: ", L₂_Q)
    # ──────────────────────────────────────────────────────────
    sheet = xf[1]
    XLSX.rename!(sheet, "circular")
    sheet["A1"] = "type w"
    sheet["B1"] = "nʷ"
    sheet["C1"] = "type φ"
    sheet["D1"] = "nᵠ"
    sheet["E1"] = "type Q"
    sheet["F1"] = "nᵛ"
    sheet["G1"] = "L₂w"
    sheet["H1"] = "L₂φ"
    sheet["I1"] = "L₂Q"
    sheet["A$row"] = "$type_w"
    sheet["B$row"] = nʷ
    sheet["C$row"] = "$type_φ"
    sheet["D$row"] = nᵠ
    sheet["E$row"] = "$type_q"
    sheet["F$row"] = nᵛ
    sheet["G$row"] = log10(L₂_w)
    sheet["H$row"] = log10(L₂_φ)
    sheet["I$row"] = log10(L₂_Q)

end
 end

gmsh.finalize()

println(to)

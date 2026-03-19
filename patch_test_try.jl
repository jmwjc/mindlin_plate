using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, L₂, L₂φ

# 新增：自然边界载荷项（与线性场一致的 M、Q）
import ApproxOperator.MindlinPlate: ∫φMdΓ, ∫wVdΓ

using TimerOutputs
import Gmsh: gmsh

# 诊断开关：默认只输出两行 L2；需要详细信息时设 PATCH_TEST_VERBOSE=1
const VERBOSE = get(ENV, "PATCH_TEST_VERBOSE", "0") in ("1", "true", "TRUE", "yes", "YES")

# ─────────────────────────────────────────────────────────────
# patch_test_try.jl（最终版输出：L₂ error of w / L₂ error of φ）
# - 线性精确解（与 patch_test.jl 一致）
# - 边界罚项采用“圆板类比映射”：
#     outer = xmax 边：w=exact, φ1=exact, φ2=exact
#     symx  = xmin 边：只强制 φ1=exact
#     symy  = ymin 边：只强制 φ2=exact
#     其余边：不施加罚项（自然边）
# - 目标：在保持线性场的前提下输出 L2 误差两行，便于对照原 patch_test.jl
# ─────────────────────────────────────────────────────────────

E = 1.0
ν = 0.3
h = 1e-1
Dᵇ = E * h^3 / 12 / (1 - ν^2)
Dˢ = (5 / 6) * E * h / (2 * (1 + ν))

# 线性精确解：选 w = a*x + b*y + c，使 φ = ∇w
# 这样 (∇w-φ)=0 ⇒ Q=0；且 φ 常量 ⇒ κ=0 ⇒ M=0；因此自然边界载荷 V=0, M₁=M₂=0
const a_w = 1.0
const b_w = -2.0
const c_w = 0.5
w(x, y, z) = a_w * x + b_w * y + c_w
w₁(x, y, z) = a_w
w₂(x, y, z) = b_w
w₁₁(x, y, z) = 0.0
w₂₂(x, y, z) = 0.0

φ₁(x, y, z) = w₁(x, y, z)
φ₂(x, y, z) = w₂(x, y, z)

# 为“纯线性场”排查：线性 φ 的二阶及更高阶导数应为 0
# 这样 κ=∇φ 为常量、弯矩 M 为 0；自然边界只剩剪力项（本例也为常量/可为 0）
φ₁₁(x, y, z) = 0.0
φ₁₂(x, y, z) = 0.0
φ₂₁(x, y, z) = 0.0
φ₂₂(x, y, z) = 0.0
φ₁₁₁(x, y, z) = 0.0
φ₁₁₂(x, y, z) = 0.0
φ₂₂₁(x, y, z) = 0.0
φ₂₂₂(x, y, z) = 0.0
φ₁₂₁(x, y, z) = 0.0
φ₁₂₂(x, y, z) = 0.0

M₁₁(x, y, z) = -Dᵇ * (φ₁₁(x, y, z) + ν * φ₂₂(x, y, z))
M₁₂(x, y, z) = -Dᵇ * (1 - ν) * 0.5 * (φ₁₂(x, y, z) + φ₂₁(x, y, z))
M₂₂(x, y, z) = -Dᵇ * (ν * φ₁₁(x, y, z) + φ₂₂(x, y, z))

M₁₁₁(x, y, z) = -Dᵇ * (φ₁₁₁(x, y, z) + ν * φ₂₂₁(x, y, z))
M₁₂₂(x, y, z) = -Dᵇ * (1 - ν) * φ₁₂₂(x, y, z)
M₁₂₁(x, y, z) = -Dᵇ * (1 - ν) * φ₁₂₁(x, y, z)
M₂₂₂(x, y, z) = -Dᵇ * (ν * φ₁₁₂(x, y, z) + φ₂₂₂(x, y, z))

Q₁(x, y, z) = Dˢ * (w₁(x, y, z) - φ₁(x, y, z))
Q₂(x, y, z) = Dˢ * (w₂(x, y, z) - φ₂(x, y, z))
Q₁₁(x, y, z) = Dˢ * (w₁₁(x, y, z) - φ₁₁(x, y, z))
Q₂₂(x, y, z) = Dˢ * (w₂₂(x, y, z) - φ₂₂(x, y, z))

# ─────────────────────────────────────────────────────────────
# 自然边界所需的边界载荷（弯矩合力 M₁,M₂ 与剪力合力 V）
# 约定与 ApproxOperator.MindlinPlate 中的 ∫φMdΓ / ∫wVdΓ 一致：
#   M₁ = n₁*M₁₁ + n₂*M₁₂
#   M₂ = n₁*M₁₂ + n₂*M₂₂
#   V  = n₁*Q₁  + n₂*Q₂
# 这里 n₁,n₂ 为边界外法向（由 GmshImport 在边界积分点提供）。
# ─────────────────────────────────────────────────────────────
# 说明：对本 patch test 选取的线性场，M=0 且 Q=0，所以任何边界自然载荷都应为 0。
# 这里保留接口函数但不再装配自然载荷项。
M₁(x, y, z, n₁, n₂) = 0.0
M₂(x, y, z, n₁, n₂) = 0.0
V(x, y, z, n₁, n₂)  = 0.0

# 说明：当前 ApproxOperator 的边界积分点未提供 :n₁/:n₂，且 prescribe! 仅支持
# f(x,y,z) / f(x,y,z,n₁) / f(x,y,z,n₁,n₂,n₃)。因此对本 patch test 的自然边(ymax)
# 直接使用“已知外法向”给出常量 M₁/M₂/V，避免依赖 ξ.n₁/ξ.n₂。
# ymax 边外法向为 (0,1)。线性场下 Q₁,Q₂,M₁₁,M₁₂,M₂₂ 要么常量要么可直接按此边求值。
const _n1_ymax = 0.0
const _n2_ymax = 1.0
M₁_ymax(x, y, z) = 0.0
M₂_ymax(x, y, z) = 0.0
V_ymax(x, y, z)  = 0.0

# 体载荷与等效弯矩载荷（域内右端项）
q(x, y, z) = -Q₁₁(x, y, z) - Q₂₂(x, y, z)
m₁(x, y, z) = M₁₁₁(x, y, z) + M₁₂₂(x, y, z) - Q₁(x, y, z)
m₂(x, y, z) = M₁₂₁(x, y, z) + M₂₂₂(x, y, z) - Q₂(x, y, z)

# 罚因子（可按需要调整；mincheck 已证明施加路径正确）
const αw_factor = 1e8
const αφ_factor = 1e8

const to = TimerOutput()

gmsh.initialize()
@timeit to "open msh file" gmsh.open("/home/jason/vscode/mindlin_plate/msh/circular_tri3_3.msh")
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()

nʷ = length(nodes)
nᵠ = length(nodes)

kʷʷ = zeros(nʷ, nʷ)
kᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)
kᵠʷ = zeros(2 * nᵠ, nʷ)
fʷ = zeros(nʷ)
fᵠ = zeros(2 * nᵠ)

@timeit to "assemble domain" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"])
    prescribe!(elements, :E => E, :ν => ν, :h => h, :q => q, :m₁ => m₁, :m₂ => m₂)
    @timeit to "calculate shape functions" set∇𝝭!(elements)

    (∫wwdΩ => elements)(kʷʷ)
    (∫φwdΩ => elements)(kᵠʷ)
    ((∫φφdΩ => elements))(kᵠᵠ)
    ((∫κκdΩ => elements))(kᵠᵠ)

    (∫wqdΩ => elements)(fʷ)
    (∫φmdΩ => elements)(fᵠ)
end

# ─────────────────────────────────────────────────────────────
# 直接按 msh 的物理组“含义”施加边界条件
# 约定：msh 中边界物理组名使用语义化命名（大小写不敏感），例如：
#   outer / xmax
#   symx  / xmin / sym-x
#   symy  / ymin / sym-y
#   natural / ymax / free
# 仍兼容旧的 Γ¹..Γ⁴：若没有语义名，则退回到几何自动分类。
# ─────────────────────────────────────────────────────────────

@timeit to "assemble boundary penalties" begin
    # 收集所有边界物理组（优先语义化命名；否则退回到所有以 Γ 开头的组）
    keylist = collect(keys(entities))
    lower_keys = Dict(lowercase(k) => k for k in keylist)

    function pick_group_names(; fallback_need_gamma::Bool=true)
        names = String[]
        # 1) 语义化关键字（大小写不敏感）
        for lk in keys(lower_keys)
            if occursin("outer", lk) || occursin("xmax", lk) || occursin("right", lk)
                push!(names, lower_keys[lk])
            elseif occursin("symx", lk) || occursin("xmin", lk) || occursin("sym-x", lk) || occursin("sym_x", lk) || occursin("left", lk)
                push!(names, lower_keys[lk])
            elseif occursin("symy", lk) || occursin("ymin", lk) || occursin("sym-y", lk) || occursin("sym_y", lk) || occursin("bottom", lk)
                push!(names, lower_keys[lk])
            elseif occursin("natural", lk) || occursin("free", lk) || occursin("ymax", lk) || occursin("top", lk)
                push!(names, lower_keys[lk])
            end
        end
        unique!(names)

        # 2) 兜底：若没有任何语义化边界组，则用 Γ* 组
        if isempty(names) && fallback_need_gamma
            names = sort([k for k in keylist if startswith(k, "Γ")])
        end
        return names
    end

    bnd_names = pick_group_names()
    if isempty(bnd_names)
        error("No boundary physical groups found (semantic or Γ*). Available keys=$(collect(keys(entities)))")
    end

    elements_Γ = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    for name in bnd_names
        elements_Γ[name] = getElements(nodes, entities[name])
        set𝝭!(elements_Γ[name])
    end

    # 映射到 4 类边界：优先按物理组名；如果是 Γ* 则用几何分类。
    Γ_outer = String[]
    Γ_symx = String[]
    Γ_symy = String[]
    Γ_natural = String[]

    used_semantic = any(!startswith(name, "Γ") for name in bnd_names)

    if used_semantic
        for name in bnd_names
            lk = lowercase(name)
            if occursin("outer", lk) || occursin("xmax", lk) || occursin("right", lk)
                push!(Γ_outer, name)
            elseif occursin("symx", lk) || occursin("xmin", lk) || occursin("sym-x", lk) || occursin("sym_x", lk) || occursin("left", lk)
                push!(Γ_symx, name)
            elseif occursin("symy", lk) || occursin("ymin", lk) || occursin("sym-y", lk) || occursin("sym_y", lk) || occursin("bottom", lk)
                push!(Γ_symy, name)
            elseif occursin("natural", lk) || occursin("free", lk) || occursin("ymax", lk) || occursin("top", lk)
                push!(Γ_natural, name)
            end
        end
    else
        side_of = Dict{String,Symbol}()
        for name in keys(elements_Γ)
            side_of[name] = classify_boundary(nodes, elements_Γ[name])
        end
        Γ_outer = [name for (name, s) in side_of if s == :xmax]
        Γ_symx = [name for (name, s) in side_of if s == :xmin]
        Γ_symy = [name for (name, s) in side_of if s == :ymin]
        Γ_natural = [name for (name, s) in side_of if s == :ymax]

        if VERBOSE
            println("[BC-map fallback] classification = ", side_of)
        end
    end

    αw = αw_factor * E
    αφ = αφ_factor * E

    if VERBOSE
        println("[BC-groups] outer=", Γ_outer, ", symx=", Γ_symx, ", symy=", Γ_symy, ", natural=", Γ_natural)
        println("[BC-param] αw=", αw, ", αφ=", αφ)
    end

    # ── outer: w + φ1, φ2
    for name in Γ_outer
        prescribe!(elements_Γ[name], :α => αw, :g => w)
        (∫αwwdΓ => elements_Γ[name])(kʷʷ, fʷ)

        prescribe!(elements_Γ[name], :α => αφ, :g₁ => φ₁, :g₂ => φ₂, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
        (∫αφφdΓ => elements_Γ[name])(kᵠᵠ, fᵠ)
    end

    # ── sym-x: only φ1
    for name in Γ_symx
        prescribe!(elements_Γ[name], :α => αφ, :g₁ => φ₁, :g₂ => φ₂, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 0.0)
        (∫αφφdΓ => elements_Γ[name])(kᵠᵠ, fᵠ)
    end

    # ── sym-y: only φ2
    for name in Γ_symy
        prescribe!(elements_Γ[name], :α => αφ, :g₁ => φ₁, :g₂ => φ₂, :n₁₁ => 0.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
        (∫αφφdΓ => elements_Γ[name])(kᵠᵠ, fᵠ)
    end

    # ── natural: assemble natural boundary loads (M, V)
    for name in Γ_natural
        prescribe!(elements_Γ[name], :M₁ => M₁_ymax, :M₂ => M₂_ymax, :V => V_ymax)
        (∫φMdΓ => elements_Γ[name])(fᵠ)
        (∫wVdΓ => elements_Γ[name])(fʷ)
    end
end

@timeit to "solve" d = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ] \ [fᵠ; fʷ]
push!(nodes, :d => d[2 * nᵠ + 1:end], :d₁ => d[1:2:2 * nᵠ], :d₂ => d[2:2:2 * nᵠ])

# 边界诊断：仅在 VERBOSE=1 时输出
if VERBOSE
    let
        # 重取边界元素（此时已有 d）
        Γall = Dict{String,Vector{ApproxOperator.AbstractElement}}()
        # 为诊断：取所有语义名 + Γ*（两者都抓）
        bnd_names_diag = sort([k for k in keys(entities) if startswith(lowercase(k), "γ") || startswith(k, "Γ") || occursin("outer", lowercase(k)) || occursin("sym", lowercase(k)) || occursin("natural", lowercase(k)) || occursin("free", lowercase(k)) || occursin("xmin", lowercase(k)) || occursin("xmax", lowercase(k)) || occursin("ymin", lowercase(k)) || occursin("ymax", lowercase(k))])
        unique!(bnd_names_diag)
        for name in bnd_names_diag
            # 防止把域 Ω 混进来
            name == "Ω" && continue
            haskey(entities, name) || continue
            Γall[name] = getElements(nodes, entities[name])
        end

        # 仍用几何分类做诊断输出（不影响施加方式）
        side_of2 = Dict{String,Symbol}()
        for name in keys(Γall)
            side_of2[name] = classify_boundary(nodes, Γall[name])
        end
        Γ_outer2 = vcat([Γall[name] for (name, s) in side_of2 if s == :xmax]...)
        Γ_symx2  = vcat([Γall[name] for (name, s) in side_of2 if s == :xmin]...)
        Γ_symy2  = vcat([Γall[name] for (name, s) in side_of2 if s == :ymin]...)
        Γ_nat2   = vcat([Γall[name] for (name, s) in side_of2 if s == :ymax]...)

        println("\n[diag-res] max residual on enforced boundaries:")
        println("  outer: max|w-w_exact|   = ", max_res_on_boundary(Γ_outer2, nd -> nd.d,  w))
        println("  outer: max|φ1-φ1_exact| = ", max_res_on_boundary(Γ_outer2, nd -> nd.d₁, φ₁))
        println("  outer: max|φ2-φ2_exact| = ", max_res_on_boundary(Γ_outer2, nd -> nd.d₂, φ₂))

        println("  symx : max|φ1-φ1_exact| = ", max_res_on_boundary(Γ_symx2,  nd -> nd.d₁, φ₁))
        println("  symy : max|φ2-φ2_exact| = ", max_res_on_boundary(Γ_symy2,  nd -> nd.d₂, φ₂))

        println("\n[observe] (not enforced) residual on other components:")
        println("  symx : max|w-w_exact|   = ", max_res_on_boundary(Γ_symx2,  nd -> nd.d,  w))
        println("  symx : max|φ2-φ2_exact| = ", max_res_on_boundary(Γ_symx2,  nd -> nd.d₂, φ₂))
        println("  symy : max|w-w_exact|   = ", max_res_on_boundary(Γ_symy2,  nd -> nd.d,  w))
        println("  symy : max|φ1-φ1_exact| = ", max_res_on_boundary(Γ_symy2,  nd -> nd.d₁, φ₁))

        println("\n[observe] natural boundary values (no penalty):")
        println("  natural: max|w-w_exact|   = ", max_res_on_boundary(Γ_nat2, nd -> nd.d,  w))
        println("  natural: max|φ1-φ1_exact| = ", max_res_on_boundary(Γ_nat2, nd -> nd.d₁, φ₁))
        println("  natural: max|φ2-φ2_exact| = ", max_res_on_boundary(Γ_nat2, nd -> nd.d₂, φ₂))
    end
end

@timeit to "calculate error" begin
    # 与 patch_test.jl 一致：用更高阶积分求 L2（这里取 10）
    @timeit to "get elements" elements_err = getElements(nodes, entities["Ω"], 10)
    prescribe!(elements_err, :E => E, :ν => ν, :h => h, :u => w, :φ₁ => φ₁, :φ₂ => φ₂)
    @timeit to "calculate shape functions" set𝝭!(elements_err)
    global L₂_w = L₂(elements_err)
    global L₂_φ = L₂φ(elements_err)
end

if VERBOSE
    # 简单观察位移/转角的数值范围
    let
        wmin = minimum(nodes.d)
        wmax = maximum(nodes.d)
        φ1min = minimum(nodes.d₁)
        φ1max = maximum(nodes.d₁)
        φ2min = minimum(nodes.d₂)
        φ2max = maximum(nodes.d₂)
        println("\n[range] w  : [", wmin, ", ", wmax, "]")
        println("[range] φ1 : [", φ1min, ", ", φ1max, "]")
        println("[range] φ2 : [", φ2min, ", ", φ2max, "]")
    end
end

gmsh.finalize()

if VERBOSE
    println(to)
end
println("L₂ error of w: ", L₂_w)
println("L₂ error of φ: ", L₂_φ)




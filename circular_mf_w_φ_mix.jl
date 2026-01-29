using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫MMdΩ, ∫∇MφdΩ, ∫MφdΓ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, ∫wVdΓ, ∫φMdΓ, L₂, L₂φ, L₂Q

using TimerOutputs, WriteVTK, XLSX
import Gmsh: gmsh

using Statistics

function cal_area_support(elms::Vector{ApproxOperator.AbstractElement})
    𝐴s = zeros(length(elms))
    for (i, elm) in enumerate(elms)
        x₁ = elm.𝓒[1].x
        y₁ = elm.𝓒[1].y
        x₂ = elm.𝓒[2].x
        y₂ = elm.𝓒[2].y
        x₃ = elm.𝓒[3].x
        y₃ = elm.𝓒[3].y
        𝐴s[i] = 0.5 * (x₁ * y₂ + x₂ * y₃ + x₃ * y₁ - x₂ * y₁ - x₃ * y₂ - x₁ * y₃)
    end
    avg𝐴 = mean(𝐴s)
    var𝐴 = var(𝐴s)
    s = (4 / 3^0.5 * avg𝐴)^0.5
    return s, var𝐴
end

E = 10.92
ν = 0.3
h = 1.0
Dᵇ = E * h^3 / 12 / (1 - ν^2)
Dˢ = 5 / 6 * E * h / (2 * (1 + ν))

R = 5.0

function w(x, y, z)
    r² = x^2 + y^2
    return (1 - r² / R^2)^2 * 11551 / 39831 * (h / 0.1)^3
end

function w₁(x, y, z)
    r² = x^2 + y^2
    return -4 * x / R^2 * (1 - r² / R^2) * 11551 / 39831 * (h / 0.1)^3
end

function w₂(x, y, z)
    r² = x^2 + y^2
    return -4 * y / R^2 * (1 - r² / R^2) * 11551 / 39831 * (h / 0.1)^3
end

function φ₁(x, y, z)
    r² = x^2 + y^2
    return -2 * x / R * (1 - r² / R^2) * 11551 / 39831 * (h / 0.1)^3
end

function φ₂(x, y, z)
    r² = x^2 + y^2
    return -2 * y / R * (1 - r² / R^2) * 11551 / 39831 * (h / 0.1)^3
end

q(x, y, z) = 1.0

Q₁(x, y, z) = Dˢ * (w₁(x, y, z) - φ₁(x, y, z))
Q₂(x, y, z) = Dˢ * (w₂(x, y, z) - φ₂(x, y, z))

M₁₁(x, y, z) = 0.0
M₁₂(x, y, z) = 0.0
M₂₂(x, y, z) = 0.0

const to = TimerOutput()

gmsh.initialize()
integrationOrder = 5
# ──────────────────────────────────────────────────────────
type_w = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_φ = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_Q = :tri3
type_M = :(PiecewisePolynomial{:Linear2D})
ndiv_φ = 16
ndiv_w = 17
ndiv = ndiv_φ

# 基于 circular.geo / circular_16.msh 的物理组命名，补齐边界别名：
#   Γ = Γᵇ ∪ Γᵉ ∪ Γˡ （全边界）
#   Γ_circ = Γᵉ       （外圆弧边界，用于固接）
function ensure_circular_boundary_aliases!(entities::Dict{String,Pair{Int,Vector{Int}}})
    # Ω 必须存在
    haskey(entities, "Ω") || error("Physical group \"Ω\" not found in msh physical names.")

    # 三段边界必须存在
    for key in ("Γᵇ", "Γᵉ", "Γˡ")
        haskey(entities, key) || error("Physical group \"$key\" not found in msh physical names.")
    end

    # 合并成 Γ
    dim = entities["Γᵇ"].first
    dim == 1 || error("Expected boundary physical groups to be dim=1, got dim=$dim for Γᵇ")
    tags = Int[]
    append!(tags, entities["Γᵇ"].second)
    append!(tags, entities["Γᵉ"].second)
    append!(tags, entities["Γˡ"].second)
    entities["Γ"] = dim => tags

    # 外圆弧边界
    entities["Γ_circ"] = entities["Γᵉ"]

    return entities
end

# NOTE: 后续如需继续修改本文件的任何逻辑（例如边界条件、对称边处理、误差指标等），将先列出：
# 1) 发现的问题/错误点；2) 产生原因；3) 建议修改方案与影响；
# 并在你确认同意后才会实际改动。

XLSX.openxlsx("xls/circular_16_tri3_17.xlsx", mode="w") do xf
    for ndiv = ndiv_w:22
        # ndiv_w = ndiv
        row = ndiv
        # ─── Deflection W ─────────────────────────────────────────
        @timeit to "open msh file" gmsh.open("msh/circular_tri3_$ndiv_w.msh")
        @timeit to "get nodes" nodes_w = get𝑿ᵢ()
        @timeit to "get entities" entities_w = getPhysicalGroups()
        ensure_circular_boundary_aliases!(entities_w)
        @timeit to "calculate support domain" begin
            elements_support = getElements(nodes_w, entities_w["Ω"], 1)
            s_w, var_A = cal_area_support(elements_support)
            γ = 1.5
            s_val = γ * s_w
            nʷ = length(nodes_w)
            push!(nodes_w, :s₁ => s_val * ones(nʷ), :s₂ => s_val * ones(nʷ), :s₃ => s_val * ones(nʷ))
        end
        xʷ = nodes_w.x
        yʷ = nodes_w.y
        zʷ = nodes_w.z
        sp_w = RegularGrid(xʷ, yʷ, zʷ, n=3, γ=5)
        # ─── Rotation Φ ───────────────────────────────────────────
        @timeit to "open msh file" gmsh.open("msh/circular_tri3_$ndiv_φ.msh")
        @timeit to "get nodes" nodes_φ = get𝑿ᵢ()
        @timeit to "get entities" entities_φ = getPhysicalGroups()
        ensure_circular_boundary_aliases!(entities_φ)
        @timeit to "calculate support domain" begin
            elements_support = getElements(nodes_φ, entities_φ["Ω"], 1)
            s_φ, var_A = cal_area_support(elements_support)
            γ = 1.5
            s_val = γ * s_φ
            nᵠ = length(nodes_φ)
            push!(nodes_φ, :s₁ => s_val * ones(nᵠ), :s₂ => s_val * ones(nᵠ), :s₃ => s_val * ones(nᵠ))
        end
        xᵠ = nodes_φ.x
        yᵠ = nodes_φ.y
        zᵠ = nodes_φ.z
        sp_φ = RegularGrid(xᵠ, yᵠ, zᵠ, n=3, γ=5)
        # ─── Shear ────────────────────────────────────────────────
        @timeit to "open msh file" gmsh.open("msh/circular_tri3_$ndiv.msh")
        @timeit to "get nodes" nodes = get𝑿ᵢ()
        @timeit to "get entities" entities = getPhysicalGroups()
        ensure_circular_boundary_aliases!(entities)

        nˢ = length(nodes)
        kᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)
        kʷʷ = zeros(nʷ, nʷ)
        kˢˢ = zeros(2 * nˢ, 2 * nˢ)
        kˢᵠ = zeros(2 * nˢ, 2 * nᵠ)
        kˢʷ = zeros(2 * nˢ, nʷ)
        kᵠʷ = zeros(2 * nᵠ, nʷ)
        fᵠ = zeros(2 * nᵠ)
        fʷ = zeros(nʷ)
        fˢ = zeros(2 * nˢ)

        @timeit to "calculate ∫QQdΩ ∫∇QwdΩ" begin
            @timeit to "get elements" elements_q = getElements(nodes, entities["Ω"], integrationOrder)
            prescribe!(elements_q, :E => E, :ν => ν, :h => h)
            @timeit to "calculate shape functions" set∇𝝭!(elements_q)

            @timeit to "get elements" elements_w = getElements(nodes_w, entities_w["Ω"], eval(type_w), integrationOrder, sp_w)
            prescribe!(elements_w, :E => E, :ν => ν, :h => h, :q => q)
            @timeit to "calculate shape functions" set𝝭!(elements_w)

            @timeit to "get elements" elements_w_Γ = getElements(nodes_w, entities_w["Γ"], eval(type_w), integrationOrder, normal=true, sp_w)
            @timeit to "calculate shape functions" set𝝭!(elements_w_Γ)

            @timeit to "get elements" elements_q_Γ = getElements(nodes, entities["Γ"], integrationOrder, normal=true)
            @timeit to "calculate shape functions" set𝝭!(elements_q_Γ)

            𝑎ˢˢ = ∫QQdΩ => elements_q
            𝑎ˢʷ = [
                ∫∇QwdΩ => (elements_q, elements_w),
                ∫QwdΓ => (elements_q_Γ, elements_w_Γ),
            ]
            𝑓ʷ = ∫wqdΩ => elements_w
            @timeit to "assemble" 𝑎ˢˢ(kˢˢ)
            @timeit to "assemble" 𝑎ˢʷ(kˢʷ)
            @timeit to "assemble" 𝑓ʷ(fʷ)
        end

        nₑ = length(elements_q)
        nᵐ = nₑ * ApproxOperator.get𝑛𝑝(eval(type_M)(𝑿ᵢ[], 𝑿ₛ[]))
        kᵐᵐ = zeros(3 * nᵐ, 3 * nᵐ)
        kᵐᵠ = zeros(3 * nᵐ, 2 * nᵠ)
        kᵐʷ = zeros(3 * nᵐ, nʷ)
        kˢᵐ = zeros(2 * nˢ, 3 * nᵐ)
        fᵐ = zeros(3 * nᵐ)

        @timeit to "calculate ∫MMdΩ ∫MφdΩ" begin
            @timeit to "get elements" elements_m = getPiecewiseElements(entities["Ω"], eval(type_M), integrationOrder)
            prescribe!(elements_m, :E => E, :ν => ν, :h => h)
            @timeit to "calculate shape functions" set∇𝝭!(elements_m)

            @timeit to "get elements" elements_φ = getElements(nodes_φ, entities_φ["Ω"], eval(type_φ), integrationOrder, sp_φ)
            prescribe!(elements_φ, :E => E, :ν => ν, :h => h)
            @timeit to "calculate shape functions" set𝝭!(elements_φ)

            @timeit to "get elements" elements_φ_Γ = getElements(nodes_φ, entities_φ["Γ"], eval(type_φ), integrationOrder, sp_φ, normal=true)
            @timeit to "calculate shape functions" set𝝭!(elements_φ_Γ)

            # =====================================================================
            # 方案A（先跑“没有纯Γ”的版本）：
            # - 不构造 piecewise 边界 elements_m_Γ
            # - 不组装边界项 ∫MφdΓ
            # 保留代码（注释）以便后续恢复。
            # =====================================================================

            # # 圆板边界由 Γᵇ、Γᵉ、Γˡ 三段组成。为避免 getPiecewiseBoundaryElements
            # # 内部对 ne/nb(=3) 的整除假设在合并后失效，这里改为分段生成再合并。
            # @timeit to "get elements" begin
            #     elements_m_Γ = ApproxOperator.AbstractElement[]
            #     append!(elements_m_Γ, getPiecewiseBoundaryElements(entities["Γᵇ"], entities["Ω"], eval(type_M), integrationOrder))
            #     append!(elements_m_Γ, getPiecewiseBoundaryElements(entities["Γᵉ"], entities["Ω"], eval(type_M), integrationOrder))
            #     append!(elements_m_Γ, getPiecewiseBoundaryElements(entities["Γˡ"], entities["Ω"], eval(type_M), integrationOrder))
            # end
            # @timeit to "calculate shape functions" set𝝭!(elements_m_Γ)

            𝑎ᵐᵐ = ∫MMdΩ => elements_m
            𝑎ᵐᵠ = [
                ∫∇MφdΩ => (elements_m, elements_φ),
                # ∫MφdΓ => (elements_m_Γ, elements_φ_Γ),
            ]
            𝑎ˢᵠ = ∫QφdΩ => (elements_q, elements_φ)
            @timeit to "assemble" 𝑎ᵐᵐ(kᵐᵐ)
            @timeit to "assemble" 𝑎ᵐᵠ(kᵐᵠ)
            @timeit to "assemble" 𝑎ˢᵠ(kˢᵠ)
        end

        # =====================================================================
        # 方案A：禁用圆弧转角罚项（依赖 elements_m_Γ）
        # =====================================================================
        # @timeit to "calculate ∫MφdΓ" begin
        #     elements_m_Γ_circ = getElements(entities["Γ_circ"], entities["Γ"], elements_m_Γ)
        #     elements_φ_Γ_circ = getElements(nodes_φ, entities_φ["Γ_circ"], eval(type_φ), integrationOrder, sp_φ, normal=true)
        #     prescribe!(elements_φ_Γ_circ, :α => 1e8 * E, :g₁ => φ₁, :g₂ => φ₂, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
        #     @timeit to "calculate shape functions" set𝝭!(elements_φ_Γ_circ)
        #     𝑎 = ∫MφdΓ => (elements_m_Γ_circ, elements_φ_Γ_circ)
        #     @timeit to "assemble" 𝑎(kᵐᵠ, fᵐ)
        # end

        @timeit to "calculate ∫QwdΓ" begin
            elements_q_Γ_circ = getElements(nodes, entities["Γ_circ"], integrationOrder, normal=true)
            elements_w_Γ_circ = getElements(nodes_w, entities_w["Γ_circ"], eval(type_w), integrationOrder, sp_w, normal=true)
            prescribe!(elements_w_Γ_circ, :α => 1e8 * E, :g => w)
            @timeit to "calculate shape functions" set𝝭!(elements_q_Γ_circ)
            @timeit to "calculate shape functions" set𝝭!(elements_w_Γ_circ)
            𝑎 = ∫QwdΓ => (elements_q_Γ_circ, elements_w_Γ_circ)
            @timeit to "assemble" 𝑎(kˢʷ, fˢ)
        end

        @timeit to "solve" d = [kᵠᵠ kᵠʷ kˢᵠ' kᵐᵠ'; kᵠʷ' kʷʷ kˢʷ' kᵐʷ'; kˢᵠ kˢʷ kˢˢ kˢᵐ; kᵐᵠ kᵐʷ kˢᵐ' kᵐᵐ] \ [fᵠ; fʷ; fˢ; fᵐ]
        push!(nodes_φ, :d₁ => d[1:2:2*nᵠ], :d₂ => d[2:2:2*nᵠ])
        push!(nodes_w, :d => d[2*nᵠ+1:2*nᵠ+nʷ])
        push!(nodes, :q₁ => d[2*nᵠ+nʷ+1:2:2*nᵠ+nʷ+2*nˢ], :q₂ => d[2*nᵠ+nʷ+2:2:2*nᵠ+nʷ+2*nˢ])
        push!(nodes, :m₁₁ => d[2*nᵠ+nʷ+2*nˢ+1:3:end], :m₂₂ => d[2*nᵠ+nʷ+2*nˢ+2:3:end], :m₁₂ => d[2*nᵠ+nʷ+2*nˢ+3:3:end])

        @timeit to "calculate error" begin
            @timeit to "get elements" elements_φ = getElements(nodes_φ, entities_φ["Ω"], eval(type_φ), 10, sp_φ)
            @timeit to "get elements" elements_w = getElements(nodes_w, entities_w["Ω"], eval(type_w), 10, sp_w)
            @timeit to "get elements" elements_q = getElements(nodes, entities["Ω"], 10)
            prescribe!(elements_φ, :E => E, :ν => ν, :h => h, :φ₁ => φ₁, :φ₂ => φ₂)
            @timeit to "calculate shape functions" set𝝭!(elements_φ)
            prescribe!(elements_w, :E => E, :ν => ν, :h => h, :u => w)
            @timeit to "calculate shape functions" set𝝭!(elements_w)
            prescribe!(elements_q, :E => E, :ν => ν, :h => h, :Q₁ => Q₁, :Q₂ => Q₂)
            @timeit to "calculate shape functions" set𝝭!(elements_q)
        end

        @timeit to "calculate error" begin
            L₂_w = L₂(elements_w)
            L₂_φ = L₂φ(elements_φ)
            L₂_Q = L₂Q(elements_q)
        end

        println(to)

        println("nʷ=$nʷ, nᵠ=$nᵠ, nˢ=$nˢ, nᵐ=$nᵐ")
        println("L₂ error of w: ", L₂_w)
        println("L₂ error of φ: ", L₂_φ)
        println("L₂ error of Q: ", L₂_Q)

        sheet = xf[1]
        XLSX.rename!(sheet, "new_sheet")
        sheet["A1"] = "type w"
        sheet["B1"] = "nʷ"
        sheet["C1"] = "type φ"
        sheet["D1"] = "nᵠ"
        sheet["E1"] = "type Q"
        sheet["F1"] = "nˢ"
        sheet["G1"] = "type M"
        sheet["H1"] = "nᵐ"
        sheet["I1"] = "L₂w"
        sheet["J1"] = "L₂φ"
        sheet["K1"] = "L₂Q"
        sheet["A$row"] = "$type_w"
        sheet["B$row"] = nʷ
        sheet["C$row"] = "$type_φ"
        sheet["D$row"] = nᵠ
        sheet["E$row"] = "$type_Q"
        sheet["F$row"] = nˢ
        sheet["G$row"] = "$type_M"
        sheet["H$row"] = nᵐ
        sheet["I$row"] = log10(L₂_w)
        sheet["J$row"] = log10(L₂_φ)
        sheet["K$row"] = log10(L₂_Q)
    end
end
gmsh.finalize()
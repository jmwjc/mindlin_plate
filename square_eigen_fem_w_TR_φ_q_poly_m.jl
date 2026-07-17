using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫Q∇wdΩ, ∫QwdΓ, ∫QφdΩ, ∫MMdΩ, ∫∇MφdΩ, ∫M∇φdΩ, ∫wwdΩ, ∫MφdΓ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, ∫wVdΓ, ∫φMdΓ, L₂w, L₂φ, L₂Q

using LinearAlgebra
using TimerOutputs, WriteVTK, XLSX 
import Gmsh: gmsh

include("cal_area_support.jl")
E = 10.92e6
ν = 0.3
h = 1e-0
Dᵇ = E/12/(1-ν^2)
Dˢ = 5/6*E/h^2/(2*(1+ν))

ndiv = 8
αʷ = 0e1*h^2
αᵠ = 0e2*h^2

const to = TimerOutput()

gmsh.initialize()

integrationOrder = 2

# type_M = :(PiecewisePolynomial{:Constant})
type_M = :(PiecewisePolynomial{:Linear2D})
# type_M = :(PiecewisePolynomial{:Quadratic2D})

@timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv.msh")
@timeit to "get nodes" nodes = get𝑿ᵢ()
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get elements" elements_φ = getElements(nodes, entities["Ω"], integrationOrder)
@timeit to "get elements" elements_Γ = getElements(nodes, entities["Γ"], integrationOrder)
elements_φ, elements_Γ = Tri3toTRTri3(elements_φ,elements_Γ)

nʷ = length(nodes)
nᵠ = length(elements_Γ)
nˢ = nᵠ
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
    @timeit to "get elements" elements_w = getElements(nodes, entities["Ω"], integrationOrder)
    prescribe!(elements_w, :E=>E, :ν=>ν, :h=>h)
    prescribe!(elements_φ, :E=>E, :ν=>ν, :h=>h)
    @timeit to "calculate shape functions" set∇𝝭!(elements_w)
    @timeit to "calculate shape functions" set∇𝝭!(elements_φ)
 
    @timeit to "get elements" elements_w_Γ = getElements(nodes, entities["Γ"], integrationOrder, normal=true)
    @timeit to "calculate shape functions" set𝝭!(elements_w_Γ)
    @timeit to "get elements" elements_φ_Γ = getElements(nodes, entities["Γ"], integrationOrder)
    elements_φ_Γ = Seg2toTRTri3(elements_φ_Γ,elements_φ)
    @timeit to "calculate shape functions" set∇𝝭!(elements_φ_Γ)

    𝑎ˢˢ = ∫QQdΩ=>elements_φ
    𝑎ˢʷ = ∫Q∇wdΩ=>(elements_φ,elements_w)
    # 𝑎ˢʷ = [
    #     ∫∇QwdΩ=>(elements_φ,elements_w),
    #     ∫QwdΓ=>(elements_φ_Γ,elements_w_Γ),
    # ]
    𝑎ˢᵠ = ∫QφdΩ=>elements_φ
    𝑎ʷʷ = ∫wwdΩ=>elements_w
    𝑓ʷ = ∫wqdΩ=>elements_w
    𝑓ᵠ = ∫φmdΩ=>elements_φ
    @timeit to "assemble" 𝑎ˢˢ(kˢˢ)
    @timeit to "assemble" 𝑎ˢʷ(kˢʷ)
    @timeit to "assemble" 𝑎ˢᵠ(kˢᵠ)
    @timeit to "assemble" 𝑎ʷʷ(kʷʷ)
end

nₑ = length(elements_w)
nᵐ = nₑ*ApproxOperator.get𝑛𝑝(eval(type_M)(𝑿ᵢ[],𝑿ₛ[]))
kᵐᵐ = zeros(3*nᵐ,3*nᵐ)
kᵐᵠ = zeros(3*nᵐ,2*nᵠ)
kᵐʷ = zeros(3*nᵐ,nʷ)
kˢᵐ = zeros(2*nˢ,3*nᵐ)
fᵐ = zeros(3*nᵐ)

@timeit to "calculate ∫MMdΩ ∫MφdΩ" begin
    @timeit to "get elements" elements_m = getPiecewiseElements(entities["Ω"], eval(type_M), integrationOrder)
    prescribe!(elements_m, :E=>E, :ν=>ν, :h=>h)
    @timeit to "calculate shape functions" set∇𝝭!(elements_m)

    @timeit to "get elements" elements_m_Γ = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], eval(type_M), integrationOrder)
    @timeit to "calculate shape functions" set𝝭!(elements_m_Γ)

    𝑎ᵐᵐ = ∫MMdΩ=>elements_m
    𝑎ᵐᵠ = ∫M∇φdΩ=>(elements_m,elements_φ)
    # 𝑎ᵐᵠ = [
    #     ∫∇MφdΩ=>(elements_m,elements_φ),
    #     ∫MφdΓ=>(elements_m_Γ,elements_φ_Γ),
    # ]
    @timeit to "assemble" 𝑎ᵐᵐ(kᵐᵐ)
    @timeit to "assemble" 𝑎ᵐᵠ(kᵐᵠ)
end

@timeit to "calculate ∫QwdΓ ∫MφdΓ" begin
    @timeit to "get elements" elements_1_w = getElements(nodes, entities["Γ¹"], integrationOrder, normal=true)
    @timeit to "get elements" elements_2_w = getElements(nodes, entities["Γ²"], integrationOrder, normal=true)
    @timeit to "get elements" elements_3_w = getElements(nodes, entities["Γ³"], integrationOrder, normal=true)
    @timeit to "get elements" elements_4_w = getElements(nodes, entities["Γ⁴"], integrationOrder, normal=true)
    elements_1_φ = Seg2toTRTri3(elements_1_w,elements_φ)
    elements_2_φ = Seg2toTRTri3(elements_2_w,elements_φ)
    elements_3_φ = Seg2toTRTri3(elements_3_w,elements_φ)
    elements_4_φ = Seg2toTRTri3(elements_4_w,elements_φ)
    @timeit to "get elements" elements_1_m = getElements(entities["Γ¹"], entities["Γ"], elements_m_Γ)
    @timeit to "get elements" elements_2_m = getElements(entities["Γ²"], entities["Γ"], elements_m_Γ)
    @timeit to "get elements" elements_3_m = getElements(entities["Γ³"], entities["Γ"], elements_m_Γ)
    @timeit to "get elements" elements_4_m = getElements(entities["Γ⁴"], entities["Γ"], elements_m_Γ)
    prescribe!(elements_1_w, :α=>αʷ*Dˢ, :g=>w)
    prescribe!(elements_2_w, :α=>αʷ*Dˢ, :g=>w)
    prescribe!(elements_3_w, :α=>αʷ*Dˢ, :g=>w)
    prescribe!(elements_4_w, :α=>αʷ*Dˢ, :g=>w)
    prescribe!(elements_1_φ, :α=>αᵠ*E, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_2_φ, :α=>αᵠ*E, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_3_φ, :α=>αᵠ*E, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_4_φ, :α=>αᵠ*E, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    @timeit to "calculate shape functions" set𝝭!(elements_1_w)
    @timeit to "calculate shape functions" set𝝭!(elements_2_w)
    @timeit to "calculate shape functions" set𝝭!(elements_3_w)
    @timeit to "calculate shape functions" set𝝭!(elements_4_w)
    @timeit to "calculate shape functions" set𝝭!(elements_1_φ)
    @timeit to "calculate shape functions" set𝝭!(elements_2_φ)
    @timeit to "calculate shape functions" set𝝭!(elements_3_φ)
    @timeit to "calculate shape functions" set𝝭!(elements_4_φ)
    𝑎ˢʷ = ∫QwdΓ => (elements_1_φ∪elements_2_φ∪elements_3_φ∪elements_4_φ,elements_1_w∪elements_2_w∪elements_3_w∪elements_4_w)
    𝑎ᵐᵠ = ∫MφdΓ => (elements_1_m ∪ elements_2_m ∪ elements_3_m ∪ elements_4_m, elements_1_φ ∪ elements_2_φ ∪ elements_3_φ ∪ elements_4_φ)
    @timeit to "assemble" 𝑎ˢʷ(kˢʷ,fˢ)
    @timeit to "assemble" 𝑎ᵐᵠ(kᵐᵠ,fᵐ)
end

# ──────────────────────────────────────────────────────────
# dᵠ = zeros(2*nᵠ)
# dʷ = zeros(nʷ)
# dˢ = zeros(2*nˢ)
# dᵐ = zeros(3*nᵐ)
# for node in nodes
#     x = node.x
#     y = node.y
#     z = node.z
#     dᵠ[2*node.𝐼-1] = φ₁(x,y,z)
#     dᵠ[2*node.𝐼]   = φ₂(x,y,z)
# end

# for node in nodes_w
#     x = node.x
#     y = node.y
#     z = node.z
#     dʷ[node.𝐼] = w(x,y,z)
# end

# for node in nodes
#     x = node.x
#     y = node.y
#     z = node.z
#     dᵠ[2*node.𝐼-1] = φ₁(x,y,z)
#     dᵠ[2*node.𝐼]   = φ₂(x,y,z)
#     dᵐ[3*node.𝐼-2] = M₁₁(x,y,z)
#     dᵐ[3*node.𝐼-1] = M₂₂(x,y,z)
#     dᵐ[3*node.𝐼]   = M₁₂(x,y,z)
# end
# println(kᵐᵠ'*dᵐ)
# println(fᵠ)
# println(kᵠᵠ*dᵠ+kᵠʷ*dʷ+kˢᵠ'*dˢ+kᵐᵠ'*dᵐ - fᵠ)
# println(kᵠʷ'*dᵠ+kʷʷ*dʷ+kˢʷ'*dˢ+kᵐʷ'*dᵐ - fʷ)
# println(kˢᵠ*dᵠ+kˢʷ*dʷ+kˢˢ*dˢ+kˢᵐ*dᵐ - fˢ)
# println(kᵐᵠ*dᵠ+kᵐʷ*dʷ+kˢᵐ'*dˢ+kᵐᵐ*dᵐ - fᵐ)
# ──────────────────────────────────────────────────────────
# @timeit to "solve" d = [kᵠᵠ kᵠʷ kˢᵠ' kᵐᵠ';kᵠʷ' kʷʷ kˢʷ' kᵐʷ';kˢᵠ kˢʷ kˢˢ kˢᵐ;kᵐᵠ kᵐʷ kˢᵐ' kᵐᵐ]\[fᵠ;fʷ;fˢ;fᵐ]
# push!(nodes,:d=>d[2*nᵠ+1:2*nᵠ+nʷ])
# ──────────────────────────────────────────────────────────
println("h = $h, Dˢ = $Dˢ, Dᵇ = $Dᵇ, nᵠ = $nᵠ, nʷ = $nʷ, nˢ = $nˢ")
k̃ᵠᵠ = - kˢᵠ'*(kˢˢ\kˢᵠ)
k̃ʷʷ = - kˢʷ'*(kˢˢ\kˢʷ)
k̃ᵠʷ = - kˢᵠ'*(kˢˢ\kˢʷ)
# ─── Eigen Test For βʷ ────────────────────────────────────
print("nʷ≤⌊[nˢ]⌋-1:         ")
n = floor(0.5*((1+8*nˢ)^0.5-3))
n_diff = 0.5*n*(n+1)-nʷ
n_diff≥0.0 ? println("✓:$n_diff") : println("×:$n_diff")
# βʷ² = eigvals(-kᵐᵠ'*(kᵐᵐ\kᵐᵠ)*(1/ndiv)^(-2))
βʷ² = eigvals(kˢʷ*(kʷʷ\kˢʷ')*(1/ndiv)^(-2))
βʷ² = real.(βʷ²)
βʷ²⁺ = βʷ²[βʷ² .≥ 1e5*eps()]
βʷ⁺ = βʷ²⁺.^0.5
nʷ⁺ = length(βʷ⁺)
# println(βʷ⁺)
βʷ⁺ = min(βʷ⁺...)
println("βʷ⁺ = $βʷ⁺, nʷ⁺ = $nʷ⁺")


print("2nˢ≤nʷ+2nᵠ-min(nʷ,n):")
n = floor(0.5*((1+8*nᵠ)^0.5-3))
# n = ceil(0.5*((1+8*nᵠ)^0.5-3))
n = 0.5*(n+2)*(n+3)
n_diff = 0.5*(nʷ+2nᵠ-min(nʷ,n))-nˢ
n_diff≥0.0 ? println("✓:$n_diff") : println("×:$n_diff")
println("nʷ = $nʷ, n = $n")

βᵞ² = eigvals([k̃ᵠᵠ k̃ᵠʷ;k̃ᵠʷ' k̃ʷʷ]/Dˢ*(1/ndiv)^(-2))
βᵞ² = real.(βᵞ²)
βᵞ²⁺ = βᵞ²[βᵞ² .≥ 1e5*eps()]
βᵞ⁺ = βᵞ²⁺.^0.5
nᵞ⁺ = length(βᵞ⁺)
βᵞ⁺ = min(βᵞ⁺...)
println("βᵞ⁺ = $βᵞ⁺, nᵞ⁺ = $nᵞ⁺")

logβʷ⁺ = log10(βʷ⁺)
logβᵞ⁺ = log10(βᵞ⁺)
println("$logβʷ⁺, $logβᵞ⁺")
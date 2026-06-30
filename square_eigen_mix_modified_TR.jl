using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QφdΩ, ∫Q∇wdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QQdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫wVdΓ, ∫αwwdΓ, ∫αφφdΓ, L₂w, L₂φ, L₂Q, H₁

using TimerOutputs, LinearAlgebra
import Gmsh: gmsh

E = 10.92e6
ν = 0.3
h = 1e-3
Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))

w(x,y,z) = 0.0
φ₁(x,y,z) = 0.0
φ₂(x,y,z) = 0.0

integrationOrder = 2
ndiv = 4

const to = TimerOutput()

gmsh.initialize()
@timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv.msh")
# @timeit to "open msh file" gmsh.open("msh/patchtest_un_tri3_$ndiv.msh")
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()
@timeit to "get elements" elements_φ = getElements(nodes, entities["Ω"], integrationOrder)
@timeit to "get elements" elements_Γ = getElements(nodes, entities["Γ"], integrationOrder)
elements_φ, elements_Γ = Tri3toTRTri3(elements_φ,elements_Γ)

nʷ = length(nodes)
nᵠ = length(elements_Γ)
nˢ = nᵠ
kʷʷ = zeros(nʷ,nʷ)
kᵠᵠ = zeros(2*nᵠ,2*nᵠ)
kˢˢ = zeros(2*nˢ,2*nˢ)
kᵠʷ = zeros(2*nᵠ,nʷ)
kˢʷ = zeros(2*nˢ,nʷ)
kˢᵠ = zeros(2*nˢ,2*nᵠ)
fʷ = zeros(nʷ)
fᵠ = zeros(2*nᵠ)
fˢ = zeros(2*nˢ)

@timeit to "calculate ∫κκdΩ" begin
    @timeit to "get elements" elements_w = getElements(nodes, entities["Ω"], integrationOrder)
    @timeit to "get elements" elements_Γ = getElements(nodes, entities["Γ"], integrationOrder)
    elements_φ_Γ = Seg2toTRTri3(elements_Γ,elements_φ)
    @timeit to "get elements" elements_w_Γ = getElements(nodes, entities["Γ"], integrationOrder)
    prescribe!(elements_w, :E=>E, :ν=>ν, :h=>h)
    prescribe!(elements_φ, :E=>E, :ν=>ν, :h=>h)
    @timeit to "calculate shape functions" set∇𝝭!(elements_w)
    @timeit to "calculate shape functions" set∇𝝭!(elements_φ)
    @timeit to "calculate shape functions" set𝝭!(elements_w_Γ)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_Γ)
    𝑎ᵠᵠ = ∫κκdΩ=>elements_φ
    𝑎ˢᵠ = ∫QφdΩ=>elements_φ
    𝑎ˢˢ = ∫QQdΩ=>elements_φ
    𝑎ˢʷ = ∫Q∇wdΩ=>(elements_φ,elements_w)
    # 𝑎ˢʷ = [
    #     ∫∇QwdΩ=>(elements_φ,elements_w),
    #     ∫QwdΓ=>(elements_φ_Γ,elements_w_Γ),
    # ]
    𝑓ʷ = ∫wqdΩ=>elements_w
    𝑓ᵠ = ∫φmdΩ=>elements_φ
    @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)
    @timeit to "assemble" 𝑎ˢᵠ(kˢᵠ)
    @timeit to "assemble" 𝑎ˢˢ(kˢˢ)
    @timeit to "assemble" 𝑎ˢʷ(kˢʷ)
end

α = 1e8
@timeit to "calculate ∫αwwdΓ" begin
    @timeit to "get elements" elements_1_w = getElements(nodes, entities["Γ¹"], integrationOrder)
    @timeit to "get elements" elements_2_w = getElements(nodes, entities["Γ²"], integrationOrder)
    @timeit to "get elements" elements_3_w = getElements(nodes, entities["Γ³"], integrationOrder)
    @timeit to "get elements" elements_4_w = getElements(nodes, entities["Γ⁴"], integrationOrder)
    elements_1_φ = Seg2toTRTri3(elements_1_w,elements_φ)
    elements_2_φ = Seg2toTRTri3(elements_2_w,elements_φ)
    elements_3_φ = Seg2toTRTri3(elements_3_w,elements_φ)
    elements_4_φ = Seg2toTRTri3(elements_4_w,elements_φ)
    prescribe!(elements_1_w, :α=>α*Dˢ, :g=>w)
    prescribe!(elements_2_w, :α=>α*Dˢ, :g=>w)
    prescribe!(elements_3_w, :α=>α*Dˢ, :g=>w)
    prescribe!(elements_4_w, :α=>α*Dˢ, :g=>w)
    prescribe!(elements_1_φ, :α=>α*Dᵇ, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_2_φ, :α=>α*Dᵇ, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_3_φ, :α=>α*Dᵇ, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_4_φ, :α=>α*Dᵇ, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    @timeit to "calculate shape functions" set𝝭!(elements_1_w)
    @timeit to "calculate shape functions" set𝝭!(elements_2_w)
    @timeit to "calculate shape functions" set𝝭!(elements_3_w)
    @timeit to "calculate shape functions" set𝝭!(elements_4_w)
    @timeit to "calculate shape functions" set𝝭!(elements_1_φ)
    @timeit to "calculate shape functions" set𝝭!(elements_2_φ)
    @timeit to "calculate shape functions" set𝝭!(elements_3_φ)
    @timeit to "calculate shape functions" set𝝭!(elements_4_φ)
    𝑎 = ∫QwdΓ => (elements_1_φ∪elements_2_φ∪elements_3_φ∪elements_4_φ,elements_1_w∪elements_2_w∪elements_3_w∪elements_4_w)
    𝑎ʷ = ∫αwwdΓ=>elements_1_w∪elements_2_w∪elements_3_w∪elements_4_w
    𝑎ᵠ = ∫αφφdΓ=>elements_1_φ∪elements_2_φ∪elements_3_φ∪elements_4_φ
    @timeit to "assemble" 𝑎(kˢʷ,fˢ)
    @timeit to "assemble" 𝑎ᵠ(kᵠᵠ,fᵠ)
end

gmsh.finalize()
println("h = $h, Dˢ = $Dˢ, Dᵇ = $Dᵇ, nᵠ = $nᵠ, nʷ = $nʷ, nˢ = $nˢ")

k̃ʷʷ = - kˢʷ'*(kˢˢ\kˢʷ)
k̃ᵠʷ = - kˢᵠ'*(kˢˢ\kˢʷ)
# ─── Eigen Test For βʷ ────────────────────────────────────
print("nʷ≤⌊[nˢ]⌋-1:         ")
n = floor(0.5*((1+8*nˢ)^0.5-3))
n_diff = 0.5*n*(n+1)-nʷ
n_diff≥0.0 ? println("✓:$n_diff") : println("×:$n_diff")
βʷ² = eigvals(k̃ʷʷ)
# βʷ² = eigvals(-kˢʷ*(kʷʷ\kˢʷ'), kˢˢ)
# βʷ² = eigvals(-kˢʷ*(k̃ʷʷ\kˢʷ'),kˢˢ)
βʷ² = real.(βʷ²)
βʷ²⁺ = βʷ²[βʷ² .≥ 1e6*eps()]
βʷ⁺ = βʷ²⁺.^0.5
nʷ⁺ = length(βʷ⁺)
# println(βʷ⁺)
βʷ⁺ = min(βʷ⁺...)
println("βʷ⁺ = $βʷ⁺, nʷ⁺ = $nʷ⁺")

# ─── Eigen Test For βᵞ ────────────────────────────────────
print("2nˢ≤nʷ+2nᵠ-min(nʷ,n):")
n = floor(0.5*((1+8*nᵠ)^0.5-3))
# n = ceil(0.5*((1+8*nᵠ)^0.5-3))
n = 0.5*(n+2)*(n+3)
n_diff = 0.5*(nʷ+2nᵠ-min(nʷ,n))-nˢ
n_diff≥0.0 ? println("✓:$n_diff") : println("×:$n_diff")
println("nʷ = $nʷ, n = $n")

# βᵞ² = eigvals(k̃ᵠʷ*k̃ᵠʷ')
# βᵞ² = eigvals(-kˢᵠ'*(kˢˢ\kˢᵠ)-k̃ᵠʷ*(k̃ʷʷ\k̃ᵠʷ'),kᵠᵠ)
βᵞ² = eigvals([-kˢᵠ'*(kˢˢ\kˢᵠ) k̃ᵠʷ;k̃ᵠʷ' k̃ʷʷ],[kᵠᵠ kᵠʷ;kᵠʷ' kʷʷ])
# βᵞ² = eigvals([-kˢᵠ'*(kˢˢ\kˢᵠ) k̃ᵠʷ;k̃ᵠʷ' k̃ʷʷ])
# βᵞ² = eigvals(-kˢᵠ'*(kˢˢ\kˢᵠ),k̃ᵠᵠ)
# βᵞ² = eigvals(k̃ᵠᵠ)
βᵞ² = real.(βᵞ²)
βᵞ²⁺ = βᵞ²[βᵞ² .≥ 1e7*eps()]
βᵞ⁺ = βᵞ²⁺.^0.5
nᵞ⁺ = length(βᵞ⁺)
βᵞ⁺ = min(βᵞ⁺...)
println("βᵞ⁺ = $βᵞ⁺, nᵞ⁺ = $nᵞ⁺")

# ─── Eigen Test For βʷ+βᵞ ──────────────────────────────────
# βʷᵞ² = eigvals(-kˢᵠ'*(kˢˢ\kˢᵠ)-k̃ᵠʷ*(k̃ʷʷ\k̃ᵠʷ'),k̃ᵠᵠ)
# βʷᵞ² = real.(βʷᵞ²)
# βʷᵞ²⁺ = βʷᵞ²[βʷᵞ² .≥ 1e5*eps()]
# βʷᵞ⁺ = βᵞ²⁺.^0.5
# nʷᵞ⁺ = length(βʷᵞ⁺)
# βʷᵞ⁺ = min(βʷᵞ⁺...)
# println("βʷᵞ⁺ = $βʷᵞ⁺, nʷᵞ⁺ = $nʷᵞ⁺")
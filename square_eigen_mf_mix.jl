using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QQdΩ, ∫wwdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, ∫wVdΓ, ∫φMdΓ, L₂w, L₂φ, L₂Q

using TimerOutputs, LinearAlgebra
import Gmsh: gmsh
include("cal_area_support.jl")

E = 10.92e6
ν = 0.3
h = 1e-0
Dᵇ = E/12/(1-ν^2)
Dˢ = 5/6*E/h^2/(2*(1+ν))

const to = TimerOutput()

# ndiv = 4, nʷ = 12, 21
# ndiv = 8, nʷ = 71, 97
# ndiv = 16, nʷ = 238, 297
# ndiv = 32, nʷ = 977, 1034, 1051, 1179
ndiv = 4
ndiv_w = ndiv
# nʷ = 1034
gmsh.initialize()
# @timeit to "open msh file" gmsh.open("msh/patchtest_tri3_irregular_$nʷ.msh")
@timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv_w.msh")
@timeit to "get nodes" nodes_w = get𝑿ᵢ()
xʷ = nodes_w.x
yʷ = nodes_w.y
zʷ = nodes_w.z
nʷ = length(nodes_w)
sp = RegularGrid(xʷ,yʷ,zʷ,n = 3,γ = 5)
@timeit to "get entities" entities = getPhysicalGroups()
elements_support = getElements(nodes_w, entities["Ω"], 1)
s, var_A = cal_area_support(elements_support)
s₁ = 1.5*s*ones(nʷ)
s₂ = 1.5*s*ones(nʷ)
s₃ = 1.5*s*ones(nʷ)
push!(nodes_w,:s₁=>s₁,:s₂=>s₂,:s₃=>s₃)
@timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv.msh")
@timeit to "get nodes" nodes = get𝑿ᵢ()
@timeit to "get entities" entities = getPhysicalGroups()

type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
nʷ = length(nodes_w)
nᵠ = length(nodes)
nˢ = length(nodes)

kʷʷ = zeros(nʷ,nʷ)
kᵠᵠ = zeros(2*nᵠ,2*nᵠ)
kˢˢ = zeros(2*nˢ,2*nˢ)
kᵠʷ = zeros(2*nᵠ,nʷ)
kˢʷ = zeros(2*nˢ,nʷ)
kˢᵠ = zeros(2*nˢ,2*nᵠ)
fʷ = zeros(nʷ)
fᵠ = zeros(2*nᵠ)
fˢ = zeros(2*nˢ)

integrationOrder = 2
@timeit to "calculate ∫κκdΩ" begin
    @timeit to "get elements" elements_w = getElements(nodes_w, entities["Ω"], type, integrationOrder, sp)
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"], integrationOrder)
    @timeit to "get elements" elements_w_Γ = getElements(nodes_w, entities["Γ"], type, integrationOrder, sp, normal=true)
    @timeit to "get elements" elements_Γ = getElements(nodes, entities["Γ"], integrationOrder, normal=true)
    prescribe!(elements, :E=>E, :ν=>ν, :h=>h)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    @timeit to "calculate shape functions" set𝝭!(elements_w)
    @timeit to "calculate shape functions" set𝝭!(elements_Γ)
    @timeit to "calculate shape functions" set𝝭!(elements_w_Γ)
    𝑎ˢᵠ = ∫QφdΩ=>elements
    𝑎ˢˢ = ∫QQdΩ=>elements
    𝑎ˢʷ = [
        ∫∇QwdΩ=>(elements,elements_w),
        ∫QwdΓ=>(elements_Γ,elements_w_Γ),
    ]
    𝑎ʷʷ = ∫wwdΩ=>elements_w
    @timeit to "assemble" 𝑎ˢˢ(kˢˢ)
    @timeit to "assemble" 𝑎ˢᵠ(kˢᵠ)
    @timeit to "assemble" 𝑎ˢʷ(kˢʷ)
    @timeit to "assemble" 𝑎ʷʷ(kʷʷ)
end

@timeit to "calculate ∫αwwdΓ ∫QwdΓ" begin
    @timeit to "get elements" elements_1 = getElements(nodes, entities["Γ¹"], integrationOrder, normal=true)
    @timeit to "get elements" elements_2 = getElements(nodes, entities["Γ²"], integrationOrder, normal=true)
    @timeit to "get elements" elements_3 = getElements(nodes, entities["Γ³"], integrationOrder, normal=true)
    @timeit to "get elements" elements_4 = getElements(nodes, entities["Γ⁴"], integrationOrder, normal=true)
    @timeit to "get elements" elements_w_1 = getElements(nodes_w, entities["Γ¹"], type, integrationOrder, sp, normal=true)
    @timeit to "get elements" elements_w_2 = getElements(nodes_w, entities["Γ²"], type, integrationOrder, sp, normal=true)
    @timeit to "get elements" elements_w_3 = getElements(nodes_w, entities["Γ³"], type, integrationOrder, sp, normal=true)
    @timeit to "get elements" elements_w_4 = getElements(nodes_w, entities["Γ⁴"], type, integrationOrder, sp, normal=true)
    prescribe!(elements_w_1, :α=>αʷ*E, :g=>0.0)
    prescribe!(elements_w_2, :α=>αʷ*E, :g=>0.0)
    prescribe!(elements_w_3, :α=>αʷ*E, :g=>0.0)
    prescribe!(elements_w_4, :α=>αʷ*E, :g=>0.0)
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    @timeit to "calculate shape functions" set𝝭!(elements_4)
    @timeit to "calculate shape functions" set𝝭!(elements_w_1)
    @timeit to "calculate shape functions" set𝝭!(elements_w_2)
    @timeit to "calculate shape functions" set𝝭!(elements_w_3)
    @timeit to "calculate shape functions" set𝝭!(elements_w_4)
    𝑎ˢ = ∫QwdΓ=>(elements_1∪elements_2∪elements_3∪elements_4,elements_w_1∪elements_w_2∪elements_w_3∪elements_w_4)
    @timeit to "assemble" 𝑎ˢ(kˢʷ,fˢ)
end

gmsh.finalize()

println("h = $h, Dˢ = $Dˢ, Dᵇ = $Dᵇ, nᵠ = $nᵠ, nʷ = $nʷ, nˢ = $nˢ")

k̃ᵠᵠ = - kˢᵠ'*(kˢˢ\kˢᵠ)
k̃ʷʷ = - kˢʷ'*(kˢˢ\kˢʷ)
k̃ᵠʷ = - kˢᵠ'*(kˢˢ\kˢʷ)
# ─── Eigen Test For βʷ ────────────────────────────────────
print("nʷ≤⌊[nˢ]⌋-1:         ")
n = floor(0.5*((1+8*nˢ)^0.5-3))
n_diff = 0.5*n*(n+1)-nʷ
n_diff≥0.0 ? println("✓:$n_diff") : println("×:$n_diff")
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
βᵞ²⁺ = βᵞ²[βᵞ² .≥ 1e4*eps()]
βᵞ⁺ = βᵞ²⁺.^0.5
nᵞ⁺ = length(βᵞ⁺)
βᵞ⁺ = min(βᵞ⁺...)
println("βᵞ⁺ = $βᵞ⁺, nᵞ⁺ = $nᵞ⁺")

logβʷ⁺ = log10(βʷ⁺)
logβᵞ⁺ = log10(βᵞ⁺)
println("$logβʷ⁺, $logβᵞ⁺")
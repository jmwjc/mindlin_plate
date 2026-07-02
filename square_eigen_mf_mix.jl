using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, ∫wVdΓ, ∫φMdΓ, L₂w, L₂φ, L₂Q

using TimerOutputs, LinearAlgebra
import Gmsh: gmsh
include("cal_area_support.jl")

E = 10.92e6
ν = 0.3
h = 1e-6
Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))

w(x,y,z) = 1/3*x^3*(x-1)^3*y^3*(y-1)^3-2*h^2/(5*(1-ν))*(y^3*(y-1)^3*x*(x-1)*(5*x^2-5*x+1)+x^3*(x-1)^3*y*(y-1)*(5*y^2-5*y+1))
w₁(x,y,z) = (x-1)^2*x^2*(2*x-1)*(y-1)^3*y^3-2*h^2/(5*(1-ν))*((20*x^3-30*x^2+12*x-1)*(y-1)^3*y^3+3*(x-1)^2*x^2*(2*x-1)*(y-1)*y*(5*y^2-5*y+1))
w₂(x,y,z) = (x-1)^3*x^3*(y-1)^2*y^2*(2*y-1)-2*h^2/(5*(1-ν))*(3*(x-1)*x*(5*x^2-5*x+1)*(y-1)^2*y^2*(2*y-1)+x^3*(x-1)^3*(20*y^3-30*y^2+12*y-1))
φ₁(x,y,z) = y^3*(y-1)^3*x^2*(x-1)^2*(2*x-1)
φ₂(x,y,z) = x^3*(x-1)^3*y^2*(y-1)^2*(2*y-1)
q(x,y,z) = E*h^3/(12*(1-ν^2))*(12*y*(y-1)*(5*x^2-5*x+1)*(2*y^2*(y-1)^2+x*(x-1)*(5*y^2-5*y+1))+12*x*(x-1)*(5*y^2-5*y+1)*(2*x^2*(x-1)^2+y*(y-1)*(5*x^2-5*x+1)))

Q₁(x,y,z) = Dˢ*(w₁(x,y,z)-φ₁(x,y,z))
Q₂(x,y,z) = Dˢ*(w₂(x,y,z)-φ₂(x,y,z))

const to = TimerOutput()

# ndiv = 4, nʷ = 12, 21
# ndiv = 8, nʷ = 71, 97
# ndiv = 16, nʷ = 238, 297
# ndiv = 32, nʷ = 977, 1034, 1051, 1179
ndiv = 32
# ndiv_w = Int(ndiv/2)
nʷ = 1034
gmsh.initialize()
@timeit to "open msh file" gmsh.open("msh/patchtest_tri3_irregular_$nʷ.msh")
@timeit to "get nodes" nodes_w = get𝑿ᵢ()
xʷ = nodes_w.x
yʷ = nodes_w.y
zʷ = nodes_w.z
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
    prescribe!(elements_w, :q=>q)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    @timeit to "calculate shape functions" set𝝭!(elements_w)
    @timeit to "calculate shape functions" set𝝭!(elements_Γ)
    @timeit to "calculate shape functions" set𝝭!(elements_w_Γ)
    𝑎ᵠᵠ = ∫κκdΩ=>elements
    𝑎ˢᵠ = ∫QφdΩ=>elements
    𝑎ˢˢ = ∫QQdΩ=>elements
    𝑎ˢʷ = [
        ∫∇QwdΩ=>(elements,elements_w),
        ∫QwdΓ=>(elements_Γ,elements_w_Γ),
    ]
    𝑓ʷ = ∫wqdΩ=>elements_w
    @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)
    @timeit to "assemble" 𝑎ˢˢ(kˢˢ)
    @timeit to "assemble" 𝑎ˢᵠ(kˢᵠ)
    @timeit to "assemble" 𝑎ˢʷ(kˢʷ)
    @timeit to "assemble" 𝑓ʷ(fʷ)
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
    prescribe!(elements_1, :α=>1e8*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_2, :α=>1e8*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_3, :α=>1e8*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_4, :α=>1e8*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_w_1,:g=>w)
    prescribe!(elements_w_2,:g=>w)
    prescribe!(elements_w_3,:g=>w)
    prescribe!(elements_w_4,:g=>w)
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    @timeit to "calculate shape functions" set𝝭!(elements_4)
    @timeit to "calculate shape functions" set𝝭!(elements_w_1)
    @timeit to "calculate shape functions" set𝝭!(elements_w_2)
    @timeit to "calculate shape functions" set𝝭!(elements_w_3)
    @timeit to "calculate shape functions" set𝝭!(elements_w_4)
    𝑎ᵠ = ∫αφφdΓ=>elements_1∪elements_2∪elements_3∪elements_4
    @timeit to "assemble" 𝑎ᵠ(kᵠᵠ,fᵠ)
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
βʷ² = eigvals(k̃ʷʷ)
βʷ² = real.(βʷ²)
βʷ²⁺ = βʷ²[βʷ² .≥ 1e6*eps()]
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

# βᵞ² = eigvals(-k̃ᵠᵠ-k̃ᵠʷ*(k̃ʷʷ\k̃ᵠʷ'),kᵠᵠ)
βᵞ² = eigvals(k̃ᵠᵠ-k̃ᵠʷ*(k̃ʷʷ\k̃ᵠʷ'))
# βᵞ² = eigvals([-kˢᵠ'*(kˢˢ\kˢᵠ) k̃ᵠʷ;k̃ᵠʷ' k̃ʷʷ],[kᵠᵠ kᵠʷ;kᵠʷ' kʷʷ])
# βᵞ² = eigvals([-kˢᵠ'*(kˢˢ\kˢᵠ) k̃ᵠʷ;k̃ᵠʷ' k̃ʷʷ])
# βᵞ² = eigvals(-kˢᵠ'*(kˢˢ\kˢᵠ),k̃ᵠᵠ)
# βᵞ² = eigvals(k̃ᵠᵠ)
βᵞ² = real.(βᵞ²)
βᵞ²⁺ = βᵞ²[βᵞ² .≥ 1e7*eps()]
βᵞ⁺ = βᵞ²⁺.^0.5
nᵞ⁺ = length(βᵞ⁺)
βᵞ⁺ = min(βᵞ⁺...)
println("βᵞ⁺ = $βᵞ⁺, nᵞ⁺ = $nᵞ⁺")

logβʷ⁺ = log10(βʷ⁺)
logβᵞ⁺ = log10(βᵞ⁺)
println("$logβʷ⁺, $logβᵞ⁺")
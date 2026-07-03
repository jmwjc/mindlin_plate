using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, ∫wVdΓ, ∫φMdΓ, L₂w, L₂φ, L₂Q, L₂γ

using TimerOutputs, WriteVTK 
import Gmsh: gmsh
include("cal_area_support.jl")

E = 10.92e6
ν = 0.3
h = 1e-3
Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))

w(x,y,z) = 1/3*x^3*(x-1)^3*y^3*(y-1)^3-2*h^2/(5*(1-ν))*(y^3*(y-1)^3*x*(x-1)*(5*x^2-5*x+1)+x^3*(x-1)^3*y*(y-1)*(5*y^2-5*y+1))
w₁(x,y,z) = (x-1)^2*x^2*(2*x-1)*(y-1)^3*y^3-2*h^2/(5*(1-ν))*((20*x^3-30*x^2+12*x-1)*(y-1)^3*y^3+3*(x-1)^2*x^2*(2*x-1)*(y-1)*y*(5*y^2-5*y+1))
w₂(x,y,z) = (x-1)^3*x^3*(y-1)^2*y^2*(2*y-1)-2*h^2/(5*(1-ν))*(3*(x-1)*x*(5*x^2-5*x+1)*(y-1)^2*y^2*(2*y-1)+x^3*(x-1)^3*(20*y^3-30*y^2+12*y-1))
φ₁(x,y,z) = y^3*(y-1)^3*x^2*(x-1)^2*(2*x-1)
φ₂(x,y,z) = x^3*(x-1)^3*y^2*(y-1)^2*(2*y-1)
q(x,y,z) = E*h^3/(12*(1-ν^2))*(12*y*(y-1)*(5*x^2-5*x+1)*(2*y^2*(y-1)^2+x*(x-1)*(5*y^2-5*y+1))+12*x*(x-1)*(5*y^2-5*y+1)*(2*x^2*(x-1)^2+y*(y-1)*(5*x^2-5*x+1)))
# γ₁(x,y,z) = w₁(x,y,z)-φ₁(x,y,z)
# γ₂(x,y,z) = w₂(x,y,z)-φ₂(x,y,z)
γ₁(x,y,z) = -2*h^2/(5*(1-ν))*((20*x^3-30*x^2+12*x-1)*(y-1)^3*y^3+3*(x-1)^2*x^2*(2*x-1)*(y-1)*y*(5*y^2-5*y+1))
γ₂(x,y,z) = -2*h^2/(5*(1-ν))*(3*(x-1)*x*(5*x^2-5*x+1)*(y-1)^2*y^2*(2*y-1)+x^3*(x-1)^3*(20*y^3-30*y^2+12*y-1))
Q₁(x,y,z) = Dˢ*γ₁(x,y,z)
Q₂(x,y,z) = Dˢ*γ₂(x,y,z)

const to = TimerOutput()

# ndiv = 4, nʷ = 12, 21
# ndiv = 8, nʷ = 71, 97
# ndiv = 16, nʷ = 238, 297
# ndiv = 32, nʷ = 977, 1034, 1051, 1179
ndiv = 32
# ndiv_w = Int(ndiv/2)
ndiv_w = ndiv
# nʷ = 12

αʷ = 1e-10
αᵠ = 1e8

gmsh.initialize()
@timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv_w.msh")
# @timeit to "open msh file" gmsh.open("msh/patchtest_tri3_irregular_$nʷ.msh")
@timeit to "get nodes" nodes_w = get𝑿ᵢ()
xʷ = nodes_w.x
yʷ = nodes_w.y
zʷ = nodes_w.z
sp = RegularGrid(xʷ,yʷ,zʷ,n = 3,γ = 5)
@timeit to "get entities" entities = getPhysicalGroups()
nʷ = length(nodes_w)
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
    prescribe!(elements_1, :α=>αᵠ*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_2, :α=>αᵠ*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_3, :α=>αᵠ*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_4, :α=>αᵠ*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_w_1, :α=>αʷ*E, :g=>w)
    prescribe!(elements_w_2, :α=>αʷ*E, :g=>w)
    prescribe!(elements_w_3, :α=>αʷ*E, :g=>w)
    prescribe!(elements_w_4, :α=>αʷ*E, :g=>w)
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
    𝑎ʷ = ∫αwwdΓ=>elements_w_1∪elements_w_2∪elements_w_3∪elements_w_4
    @timeit to "assemble" 𝑎ʷ(kʷʷ,fʷ)
end

# dᵠ = zeros(2*nᵠ)
# dˢ = zeros(2*nˢ)
# dʷ = zeros(nʷ)
# for node in nodes
#     x = node.x
#     y = node.y
#     z = node.z
#     dᵠ[2*node.𝐼-1] = φ₁(x,y,z)
#     dᵠ[2*node.𝐼]   = φ₂(x,y,z)
#     dˢ[2*node.𝐼-1] = Q₁(x,y,z)
#     dˢ[2*node.𝐼]   = Q₂(x,y,z)
#     dʷ[node.𝐼] = w(x,y,z)
# end
# println(kˢˢ*dˢ)
# println(kˢʷ*dʷ)
# println(kˢˢ*dˢ + kˢʷ*dʷ)
# println(kˢʷ*ones(nʷ).-fˢ)
# println(kᵠᵠ*dᵠ + kˢᵠ'*dˢ - fᵠ)
# println(kˢˢ*dˢ + kˢᵠ*dᵠ + kˢʷ*dʷ - fˢ)
# println(kˢʷ'*dˢ + kʷʷ*dʷ - fʷ)
# println(kˢᵠ*dᵠ)
# println(kˢʷ*dʷ)
# println(kˢˢ*dˢ)
# println(kˢˢ*dˢ + kˢʷ*dʷ)

@timeit to "solve" d = [kᵠᵠ kᵠʷ kˢᵠ';kᵠʷ' kʷʷ kˢʷ';kˢᵠ kˢʷ kˢˢ]\[fᵠ;fʷ;fˢ]
push!(nodes,:d₁=>d[1:2:2*nᵠ], :d₂=>d[2:2:2*nᵠ], :q₁=>d[2*nᵠ+nʷ+1:2:end], :q₂=>d[2*nᵠ+nʷ+2:2:end])
push!(nodes_w,:d=>d[2*nᵠ+1:2*nᵠ+nʷ])

@timeit to "calculate error" begin
    @timeit to "get elements" elements_w = getElements(nodes_w, entities["Ω"], type, 10, sp)
    prescribe!(elements_w, :E=>E, :ν=>ν, :h=>h, :w=>w, :γ₁=>γ₁, :γ₂=>γ₂)
    @timeit to "calculate shape functions" set∇𝝭!(elements_w)
    L₂_w = L₂w(elements_w)
    @timeit to "get elements" elements_φ = getElements(nodes, entities["Ω"], 10)
    prescribe!(elements_φ, :E=>E, :ν=>ν, :h=>h, :φ₁=>φ₁, :φ₂=>φ₂, :Q₁=>Q₁, :Q₂=>Q₂)
    @timeit to "calculate shape functions" set𝝭!(elements_φ)
    L₂_φ = L₂φ(elements_φ)
    L₂_Q = L₂Q(elements_φ)
end

gmsh.finalize()

println("h = $h, Dˢ = $Dˢ, Dᵇ = $Dᵇ, nᵠ = $nᵠ, nʷ = $nʷ, nˢ = $nˢ")
print("nˢ≤ᵠ:         ")
n_diff = nᵠ-nˢ
n_diff≥0.0 ? println("✓:$n_diff") : println("×:$n_diff")
print("nʷ≤⌊[nˢ]⌋-1:  ")
n = floor(0.5*((1+8*nˢ)^0.5-3))
n_diff = 0.5*n*(n+1)-nʷ
n_diff≥0.0 ? println("✓:$n_diff") : println("×:$n_diff")
println(0.5*n*(n+1))

logL₂w = log10(L₂_w)
logL₂φ = log10(L₂_φ)
logL₂Q = log10(L₂_Q)
println("$logL₂w, $logL₂φ, $logL₂Q")





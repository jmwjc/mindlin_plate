using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫MMdΩ, ∫∇MφdΩ, ∫MφdΓ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, ∫wVdΓ, ∫φMdΓ, L₂w, L₂φ, L₂Q

using LinearAlgebra
using TimerOutputs, WriteVTK, XLSX 
import Gmsh: gmsh

include("cal_area_support.jl")
E = 10.92e6
ν = 0.3
h = 1e-0
Dᵇ = E/12/(1-ν^2)
Dˢ = 5/6*E/h^2/(2*(1+ν))

ndiv_φ = 16
αʷ = 0e1*h^2
αᵠ = 0e1*h^2

w(x,y,z) = 1/3*x^3*(x-1)^3*y^3*(y-1)^3-2*h^2/(5*(1-ν))*(y^3*(y-1)^3*x*(x-1)*(5*x^2-5*x+1)+x^3*(x-1)^3*y*(y-1)*(5*y^2-5*y+1))
w₁(x,y,z) = (x-1)^2*x^2*(2*x-1)*(y-1)^3*y^3-2*h^2/(5*(1-ν))*((20*x^3-30*x^2+12*x-1)*(y-1)^3*y^3+3*(x-1)^2*x^2*(2*x-1)*(y-1)*y*(5*y^2-5*y+1))
w₂(x,y,z) = (x-1)^3*x^3*(y-1)^2*y^2*(2*y-1)-2*h^2/(5*(1-ν))*(3*(x-1)*x*(5*x^2-5*x+1)*(y-1)^2*y^2*(2*y-1)+x^3*(x-1)^3*(20*y^3-30*y^2+12*y-1))
φ₁(x,y,z) = y^3*(y-1)^3*x^2*(x-1)^2*(2*x-1)
φ₂(x,y,z) = x^3*(x-1)^3*y^2*(y-1)^2*(2*y-1)
φ₁₁(x,y,z) = y^3*(y-1)^3 * (2*x*(x-1)*(2*x-1)^2 + 2*x^2*(x-1)^2)
φ₁₂(x,y,z) = 3*y^2*(y-1)^2*(2*y-1) * (x^2*(x-1)^2*(2*x-1))
φ₂₁(x,y,z) = 3*x^2*(x-1)^2*(2*x-1) * (y^2*(y-1)^2*(2*y-1))
φ₂₂(x,y,z) = x^3*(x-1)^3 * (2*y*(y-1)*(2*y-1)^2 + 2*y^2*(y-1)^2)
q(x,y,z) = E/(12*(1-ν^2))*(12*y*(y-1)*(5*x^2-5*x+1)*(2*y^2*(y-1)^2+x*(x-1)*(5*y^2-5*y+1))+12*x*(x-1)*(5*y^2-5*y+1)*(2*x^2*(x-1)^2+y*(y-1)*(5*x^2-5*x+1)))

Q₁(x,y,z) = Dˢ*(w₁(x,y,z)-φ₁(x,y,z))
Q₂(x,y,z) = Dˢ*(w₂(x,y,z)-φ₂(x,y,z))

M₁₁(x,y,z)= -Dᵇ*(φ₁₁(x,y,z)+ν*φ₂₂(x,y,z))
M₁₂(x,y,z)= -Dᵇ*(1-ν)*0.5*(φ₁₂(x,y,z)+φ₂₁(x,y,z))
M₂₂(x,y,z)= -Dᵇ*(ν*φ₁₁(x,y,z)+φ₂₂(x,y,z))

# w(x,y,z) = - (x^4*y - x*y^4 + x^3*y^2 - x^2*y^3)/6/Dᵇ/h^3 + (x^3 + 3*x^2*y - 3*x*y^2 - y^3)/3/Dˢ/h^3
# w₁(x,y,z) = - (4*x^3*y - y^4 + 3*x^2*y^2 - 2*x*y^3)/6/Dᵇ/h^3 + (x^2 + 2*x*y - y^2)/Dˢ/h^3
# w₂(x,y,z) = - (x^4 - 4*x*y^3 + 2*x^3*y - 3*x^2*y^2)/6/Dᵇ/h^3 + (x^2 - 2*x*y - y^2)/Dˢ/h^3
# w₁₁(x,y,z) = - (6*x^2*y + 3*x*y^2 - y^3)/3/Dᵇ/h^3 + 2*(x + y)/Dˢ/h^3
# w₂₂(x,y,z) = - (- 6*x*y^2 + x^3 - 3*x^2*y)/3/Dᵇ/h^3 + 2*(- x - y)/Dˢ/h^3
# φ₁(x,y,z) = - (4*x^3*y - y^4 + 3*x^2*y^2 - 2*x*y^3)/6/Dᵇ/h^3
# φ₂(x,y,z) = - (x^4 - 4*x*y^3 + 2*x^3*y - 3*x^2*y^2)/6/Dᵇ/h^3
# φ₁₁(x,y,z) = - (6*x^2*y + 3*x*y^2 - y^3)/3/Dᵇ/h^3
# φ₁₂(x,y,z) = - (2*x^3 - 2*y^3 + 3*x^2*y - 3*x*y^2)/3/Dᵇ/h^3
# φ₂₁(x,y,z) = - (2*x^3 - 2*y^3 + 3*x^2*y - 3*x*y^2)/3/Dᵇ/h^3
# φ₂₂(x,y,z) = - (- 6*x*y^2 + x^3 - 3*x^2*y)/3/Dᵇ/h^3
# φ₁₁₁(x,y,z) = - (4*x*y + y^2)/Dᵇ/h^3
# φ₁₁₂(x,y,z) = - (2*x^2 + 2*x*y - y^2)/Dᵇ/h^3
# φ₂₂₁(x,y,z) = - (- 2*y^2 + x^2 - 2*x*y)/Dᵇ/h^3
# φ₂₂₂(x,y,z) = - (- 4*x*y - x^2)/Dᵇ/h^3
# φ₁₂₁(x,y,z) = - (2*x^2 + 2*x*y - y^2)/Dᵇ/h^3
# φ₁₂₂(x,y,z) = - (- 2*y^2 + x^2 - 2*x*y)/Dᵇ/h^3

# M₁₁(x,y,z)= -Dᵇ*(φ₁₁(x,y,z)+ν*φ₂₂(x,y,z))
# M₁₂(x,y,z)= -Dᵇ*(1-ν)*0.5*(φ₁₂(x,y,z)+φ₂₁(x,y,z))
# M₂₂(x,y,z)= -Dᵇ*(ν*φ₁₁(x,y,z)+φ₂₂(x,y,z))
# M₁₁₁(x,y,z)= -Dᵇ*(φ₁₁₁(x,y,z)+ν*φ₂₂₁(x,y,z))
# M₁₂₂(x,y,z)= -Dᵇ*(1-ν)*φ₁₂₂(x,y,z)
# M₁₂₁(x,y,z)= -Dᵇ*(1-ν)*φ₁₂₁(x,y,z)
# M₂₂₂(x,y,z)= -Dᵇ*(ν*φ₁₁₂(x,y,z)+φ₂₂₂(x,y,z))

# γ₁(x,y,z) = w₁(x,y,z) - φ₁(x,y,z)
# γ₂(x,y,z) = w₂(x,y,z) - φ₂(x,y,z)
# Q₁(x,y,z) = Dˢ*γ₁(x,y,z)
# Q₂(x,y,z) = Dˢ*γ₂(x,y,z)
# # Q₁(x,y,z) = (x^2 + 2*x*y - y^2)/h^3
# # Q₂(x,y,z) = (x^2 - 2*x*y - y^2)/h^3
# Q₁₁(x,y,z) =  2*(x + y)
# Q₂₂(x,y,z) = -2*(x + y)
# # q(x,y,z) = 0.0
# q(x,y,z)=-Q₁₁(x,y,z)-Q₂₂(x,y,z)
# m₁(x,y,z) = M₁₁₁(x,y,z)+M₁₂₂(x,y,z) - Q₁(x,y,z)
# m₂(x,y,z) = M₁₂₁(x,y,z)+M₂₂₂(x,y,z) - Q₂(x,y,z)


const to = TimerOutput()

gmsh.initialize()
# @timeit to "open msh file" gmsh.open("msh/patchtest_3.msh")
# @timeit to "get nodes" nodes_s = get𝑿ᵢ()

integrationOrder = 2
# ──────────────────────────────────────────────────────────
type_w = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_φ = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_Q = :tri3
# type_M = :(PiecewisePolynomial{:Constant})
type_M = :(PiecewisePolynomial{:Linear2D})
# type_M = :(PiecewisePolynomial{:Quadratic2D})

sʷ = 1.5
sᵠ = 1.5
XLSX.openxlsx("xls/square_$(ndiv_φ).xlsx", mode="w") do xf
row = 1
for ndiv_w = 2:24
    for ndiv = 2:24
        row += 1

# ─── Deflection W ─────────────────────────────────────────
@timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv_w.msh")
# @timeit to "open msh file" gmsh.open("msh/patchtest_tri3_irregular_$nʷ.msh")
@timeit to "get nodes" nodes_w = get𝑿ᵢ()
@timeit to "get entities" entities_w = getPhysicalGroups()
xʷ = nodes_w.x
yʷ = nodes_w.y
zʷ = nodes_w.z
sp_w = RegularGrid(xʷ,yʷ,zʷ,n = 3,γ = 5)
nʷ = length(nodes_w)
# s = 1/ndiv_w
# s₁ = sʷ * s * ones(nʷ)
# s₂ = sʷ * s * ones(nʷ)
# s₃ = sʷ * s * ones(nʷ)
# push!(nodes_w,:s₁=>s₁,:s₂=>s₂,:s₃=>s₃)
elements_support = getElements(nodes_w, entities_w["Ω"], 1)
s_w, var_A = cal_area_support(elements_support)
s = sʷ*s_w*ones(nʷ)
push!(nodes_w, :s₁=>s, :s₂=>s, :s₃=>s)
# ─── Rotation Φ ───────────────────────────────────────────
@timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv_φ.msh")
@timeit to "get nodes" nodes_φ = get𝑿ᵢ()
@timeit to "get entities" entities_φ = getPhysicalGroups()
xᵠ = nodes_φ.x
yᵠ = nodes_φ.y
zᵠ = nodes_φ.z
sp_φ = RegularGrid(xᵠ,yᵠ,zᵠ,n = 3,γ = 5)
nᵠ = length(nodes_φ)
# s = 1/ndiv_φ
# s₁ = sᵠ * s * ones(nᵠ)
# s₂ = sᵠ * s * ones(nᵠ)
# s₃ = sᵠ * s * ones(nᵠ)
# push!(nodes_φ,:s₁=>s₁,:s₂=>s₂,:s₃=>s₃)
elements_support = getElements(nodes_φ, entities_φ["Ω"], 1)
s_φ, var_A = cal_area_support(elements_support)
s = sᵠ*s_φ*ones(nᵠ)
push!(nodes_φ, :s₁=>s, :s₂=>s, :s₃=>s)
# ─── Shear ────────────────────────────────────────────────
@timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv.msh")
# @timeit to "open msh file" gmsh.open("msh/patchtest_high_un_tri3_$ndiv.msh")
@timeit to "get nodes" nodes = get𝑿ᵢ()
@timeit to "get entities" entities = getPhysicalGroups()

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
kᵅᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)
kᵅʷʷ = zeros(nʷ, nʷ)
fᵅᵠ = zeros(2 * nᵠ)
fᵅʷ = zeros(nʷ)

@timeit to "calculate ∫QQdΩ ∫∇QwdΩ" begin
    @timeit to "get elements" elements_q = getElements(nodes, entities["Ω"], integrationOrder)
    prescribe!(elements_q, :E=>E, :ν=>ν, :h=>h)
    @timeit to "calculate shape functions" set∇𝝭!(elements_q)

    @timeit to "get elements" elements_w = getElements(nodes_w, entities["Ω"], eval(type_w), integrationOrder, sp_w)
    prescribe!(elements_w, :E=>E, :ν=>ν, :h=>h, :q=>q)
    @timeit to "calculate shape functions" set𝝭!(elements_w)

    @timeit to "get elements" elements_w_Γ = getElements(nodes_w, entities["Γ"], eval(type_w), integrationOrder, sp_w, normal=true)
    @timeit to "calculate shape functions" set𝝭!(elements_w_Γ)

    @timeit to "get elements" elements_q_Γ = getElements(nodes, entities["Γ"], integrationOrder, normal=true)
    @timeit to "calculate shape functions" set𝝭!(elements_q_Γ)
 
    𝑎ˢˢ = ∫QQdΩ=>elements_q
    𝑎ˢʷ = [
        ∫∇QwdΩ=>(elements_q,elements_w),
        ∫QwdΓ=>(elements_q_Γ,elements_w_Γ),
    ]
    𝑓ʷ = ∫wqdΩ=>elements_w
    @timeit to "assemble" 𝑎ˢˢ(kˢˢ)
    @timeit to "assemble" 𝑎ˢʷ(kˢʷ)
    @timeit to "assemble" 𝑓ʷ(fʷ)
end

nₑ = length(elements_q)
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

    @timeit to "get elements" elements_φ = getElements(nodes_φ, entities["Ω"], eval(type_φ), integrationOrder, sp_φ)
    prescribe!(elements_φ, :E=>E, :ν=>ν, :h=>h)
    @timeit to "calculate shape functions" set𝝭!(elements_φ)

    @timeit to "get elements" elements_φ_Γ = getElements(nodes_φ, entities["Γ"], eval(type_φ), integrationOrder, sp_φ, normal=true)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_Γ)

    @timeit to "get elements" elements_m_Γ = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], eval(type_M), integrationOrder)
    @timeit to "calculate shape functions" set𝝭!(elements_m_Γ)

    𝑎ᵐᵐ = ∫MMdΩ=>elements_m
    𝑎ᵐᵠ = [
        ∫∇MφdΩ=>(elements_m,elements_φ),
        ∫MφdΓ=>(elements_m_Γ,elements_φ_Γ),
    ]
    𝑎ˢᵠ = ∫QφdΩ=>(elements_q,elements_φ)
    # 𝑓ᵠ = ∫φmdΩ=>elements_φ
    @timeit to "assemble" 𝑎ᵐᵐ(kᵐᵐ)
    @timeit to "assemble" 𝑎ᵐᵠ(kᵐᵠ)
    @timeit to "assemble" 𝑎ˢᵠ(kˢᵠ)
    # @timeit to "assemble" 𝑓ᵠ(fᵠ)
end

@timeit to "calculate ∫QwdΓ" begin
    @timeit to "get elements" elements_q_1 = getElements(nodes, entities["Γ¹"], integrationOrder, normal=true)
    @timeit to "get elements" elements_q_2 = getElements(nodes, entities["Γ²"], integrationOrder, normal=true)
    @timeit to "get elements" elements_q_3 = getElements(nodes, entities["Γ³"], integrationOrder, normal=true)
    @timeit to "get elements" elements_q_4 = getElements(nodes, entities["Γ⁴"], integrationOrder, normal=true)
    @timeit to "get elements" elements_w_1 = getElements(nodes_w, entities["Γ¹"], eval(type_w), integrationOrder, sp_w, normal=true)
    @timeit to "get elements" elements_w_2 = getElements(nodes_w, entities["Γ²"], eval(type_w), integrationOrder, sp_w, normal=true)
    @timeit to "get elements" elements_w_3 = getElements(nodes_w, entities["Γ³"], eval(type_w), integrationOrder, sp_w, normal=true)
    @timeit to "get elements" elements_w_4 = getElements(nodes_w, entities["Γ⁴"], eval(type_w), integrationOrder, sp_w, normal=true)
    prescribe!(elements_w_1, :α=>αʷ*Dˢ, :g=>w)
    prescribe!(elements_w_2, :α=>αʷ*Dˢ, :g=>w)
    prescribe!(elements_w_3, :α=>αʷ*Dˢ, :g=>w)
    prescribe!(elements_w_4, :α=>αʷ*Dˢ, :g=>w)
    @timeit to "calculate shape functions" set𝝭!(elements_q_1)
    @timeit to "calculate shape functions" set𝝭!(elements_q_2)
    @timeit to "calculate shape functions" set𝝭!(elements_q_3)
    @timeit to "calculate shape functions" set𝝭!(elements_q_4)
    @timeit to "calculate shape functions" set𝝭!(elements_w_1)
    @timeit to "calculate shape functions" set𝝭!(elements_w_2)
    @timeit to "calculate shape functions" set𝝭!(elements_w_3)
    @timeit to "calculate shape functions" set𝝭!(elements_w_4)
        𝑎 = ∫QwdΓ => (elements_q_1 ∪ elements_q_2 ∪ elements_q_3 ∪ elements_q_4, elements_w_1 ∪ elements_w_2 ∪ elements_w_3 ∪ elements_w_4)
    𝑎ʷ = ∫αwwdΓ => elements_w_1 ∪ elements_w_2 ∪ elements_w_3 ∪ elements_w_4
    @timeit to "assemble" 𝑎(kˢʷ,fˢ)
    @timeit to "assemble" 𝑎ʷ(kʷʷ, fʷ)
end

@timeit to "calculate ∫MφdΓ" begin
    @timeit to "get elements" elements_m_1 = getElements(entities["Γ¹"], entities["Γ"], elements_m_Γ)
    @timeit to "get elements" elements_m_2 = getElements(entities["Γ²"], entities["Γ"], elements_m_Γ)
    @timeit to "get elements" elements_m_3 = getElements(entities["Γ³"], entities["Γ"], elements_m_Γ)
    @timeit to "get elements" elements_m_4 = getElements(entities["Γ⁴"], entities["Γ"], elements_m_Γ)
    @timeit to "get elements" elements_φ_1 = getElements(nodes_φ, entities["Γ¹"], eval(type_φ), integrationOrder, sp_φ, normal=true)
    @timeit to "get elements" elements_φ_2 = getElements(nodes_φ, entities["Γ²"], eval(type_φ), integrationOrder, sp_φ, normal=true)
    @timeit to "get elements" elements_φ_3 = getElements(nodes_φ, entities["Γ³"], eval(type_φ), integrationOrder, sp_φ, normal=true)
    @timeit to "get elements" elements_φ_4 = getElements(nodes_φ, entities["Γ⁴"], eval(type_φ), integrationOrder, sp_φ, normal=true)
    prescribe!(elements_φ_1, :α=>αᵠ*E, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_φ_2, :α=>αᵠ*E, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_φ_3, :α=>αᵠ*E, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_φ_4, :α=>αᵠ*E, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_1)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_2)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_3)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_4)
    𝑎 = ∫MφdΓ => (elements_m_1 ∪ elements_m_2 ∪ elements_m_3 ∪ elements_m_4, elements_φ_1 ∪ elements_φ_2 ∪ elements_φ_3 ∪ elements_φ_4)
    𝑎ᵅ = ∫αφφdΓ => elements_φ_1 ∪ elements_φ_2 ∪ elements_φ_3 ∪ elements_φ_4
    @timeit to "assemble" 𝑎(kᵐᵠ,fᵐ)
    @timeit to "assemble" 𝑎ᵅ(kᵠᵠ, fᵠ)
end

# ──────────────────────────────────────────────────────────
# dᵠ = zeros(2*nᵠ)
# dʷ = zeros(nʷ)
# dˢ = zeros(2*nˢ)
# dᵐ = zeros(3*nᵐ)
# for node in nodes_φ
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
logL₂w = 0.0
logL₂φ = 0.0
logL₂Q = 0.0
    try
        d = [kᵠᵠ kᵠʷ kˢᵠ' kᵐᵠ';kᵠʷ' kʷʷ kˢʷ' kᵐʷ';kˢᵠ kˢʷ kˢˢ kˢᵐ;kᵐᵠ kᵐʷ kˢᵐ' kᵐᵐ]\[fᵠ;fʷ;fˢ;fᵐ]
        push!(nodes_φ,:d₁=>d[1:2:2*nᵠ], :d₂=>d[2:2:2*nᵠ])
        push!(nodes_w,:d=>d[2*nᵠ+1:2*nᵠ+nʷ])
        push!(nodes,:q₁=>d[2*nᵠ+nʷ+1:2:2*nᵠ+nʷ+2*nˢ], :q₂=>d[2*nᵠ+nʷ+2:2:2*nᵠ+nʷ+2*nˢ])
        push!(nodes,:m₁₁=>d[2*nᵠ+nʷ+2*nˢ+1:3:end],:m₂₂=>d[2*nᵠ+nʷ+2*nˢ+2:3:end],:m₁₂=>d[2*nᵠ+nʷ+2*nˢ+3:3:end])
        # ──────────────────────────────────────────────────────────
        @timeit to "calculate error" begin
            @timeit to "get elements" elements_φ = getElements(nodes_φ, entities["Ω"], eval(type_φ), 10, sp_φ)
            @timeit to "get elements" elements_w = getElements(nodes_w, entities["Ω"], eval(type_w), 10, sp_w)
            @timeit to "get elements" elements_q = getElements(nodes, entities["Ω"], 10)
            # @timeit to "get elements" elements_m = getElements(nodes, entities["Ω"], 10)
            prescribe!(elements_φ, :E=>E, :ν=>ν, :h=>h, :φ₁=>φ₁, :φ₂=>φ₂)
            @timeit to "calculate shape functions" set𝝭!(elements_φ)
            prescribe!(elements_w, :E=>E, :ν=>ν, :h=>h, :w=>w)
            @timeit to "calculate shape functions" set𝝭!(elements_w)
            prescribe!(elements_q, :E=>E, :ν=>ν, :h=>h, :Q₁=>Q₁, :Q₂=>Q₂)
            @timeit to "calculate shape functions" set𝝭!(elements_q)
            # prescribe!(elements_m, :E=>E, :ν=>ν, :h=>h, :M₁₁=>Q₁, :Q₂=>Q₂)
            # @timeit to "calculate shape functions" set𝝭!(elements_q)
        end

        @timeit to "calculate error" begin
            L₂_w = L₂w(elements_w)
            L₂_φ = L₂φ(elements_φ)
            L₂_Q = L₂Q(elements_q)
        end
        logL₂w = log10(L₂_w)
        logL₂φ = log10(L₂_φ)
        logL₂Q = log10(L₂_Q)
    catch e
        if isa(e, SingularException) || isa(e, LAPACKException)
            logL₂w = nothing
            logL₂φ = nothing
            logL₂Q = nothing
        else
            rethrow(e) 
        end
    end
# ──────────────────────────────────────────────────────────

            sheet = xf[1]
            XLSX.rename!(sheet, "new_sheet")
            sheet["A1"] = "nʷ"
            sheet["B1"] = "nˢ"
            sheet["C1"] = "L₂w"
            sheet["D1"] = "L₂φ"
            sheet["E1"] = "L₂Q"
            sheet["A$row"] = nʷ
            sheet["B$row"] = nˢ
            sheet["C$row"] = logL₂w
            sheet["D$row"] = logL₂φ
            sheet["E$row"] = logL₂Q
        end
    end
end
gmsh.finalize()



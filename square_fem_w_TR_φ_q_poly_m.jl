using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫Q∇wdΩ, ∫QwdΓ, ∫QφdΩ, ∫MMdΩ, ∫∇MφdΩ, ∫M∇φdΩ, ∫MφdΓ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, ∫wVdΓ, ∫φMdΓ, L₂w, L₂φ, L₂Q

using LinearAlgebra
using TimerOutputs, WriteVTK, XLSX 
import Gmsh: gmsh

include("cal_area_support.jl")
E = 10.92e6
ν = 0.3
h = 1e-3
Dᵇ = E/12/(1-ν^2)
Dˢ = 5/6*E/h^2/(2*(1+ν))

ndiv = 4
αʷ = 0e1*h^2
αᵠ = 0e2*h^2

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

# M₁₁(x,y,z)= -Dᵇ*(φ₁₁(x,y,z)+ν*φ₂₂(x,y,z))
# M₁₂(x,y,z)= -Dᵇ*(1-ν)*0.5*(φ₁₂(x,y,z)+φ₂₁(x,y,z))
# M₂₂(x,y,z)= -Dᵇ*(ν*φ₁₁(x,y,z)+φ₂₂(x,y,z))

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

# w(x,y,z) = 1.0+x+y
# w₁(x,y,z) = 1.0
# w₂(x,y,z) = 1.0
# w₁₁(x,y,z) = 0.0
# w₂₂(x,y,z) = 0.0

# φ₁(x,y,z) = 1.0+x+y
# φ₂(x,y,z) = 1.0+x+y
# φ₁₁(x,y,z)  = 1.0
# φ₁₂(x,y,z)  = 1.0
# φ₂₁(x,y,z)  = 1.0
# φ₂₂(x,y,z)  = 1.0

# φ₁₁₁(x,y,z)  = 0.0
# φ₁₁₂(x,y,z)  = 0.0
# φ₂₂₁(x,y,z)  = 0.0
# φ₂₂₂(x,y,z)  = 0.0
# φ₁₂₁(x,y,z)  = 0.0
# φ₁₂₂(x,y,z)  = 0.0

# M₁₁(x,y,z)= -Dᵇ*(φ₁₁(x,y,z)+ν*φ₂₂(x,y,z))
# M₁₂(x,y,z)= -Dᵇ*(1-ν)*0.5*(φ₁₂(x,y,z)+φ₂₁(x,y,z))
# M₂₂(x,y,z)= -Dᵇ*(ν*φ₁₁(x,y,z)+φ₂₂(x,y,z))
# M₁₁₁(x,y,z)= -Dᵇ*(φ₁₁₁(x,y,z)+ν*φ₂₂₁(x,y,z))
# M₁₂₂(x,y,z)= -Dᵇ*(1-ν)*φ₁₂₂(x,y,z)
# M₁₂₁(x,y,z)= -Dᵇ*(1-ν)*φ₁₂₁(x,y,z)
# M₂₂₂(x,y,z)= -Dᵇ*(ν*φ₁₁₂(x,y,z)+φ₂₂₂(x,y,z))

# Q₁(x,y,z) = Dˢ*(w₁(x,y,z)-φ₁(x,y,z))
# Q₂(x,y,z) = Dˢ*(w₂(x,y,z)-φ₂(x,y,z))
# Q₁₁(x,y,z) = Dˢ*(w₁₁(x,y,z)-φ₁₁(x,y,z))
# Q₂₂(x,y,z) = Dˢ*(w₂₂(x,y,z)-φ₂₂(x,y,z))
# q(x,y,z)=-Q₁₁(x,y,z)-Q₂₂(x,y,z)
# m₁(x,y,z) = M₁₁₁(x,y,z)+M₁₂₂(x,y,z) - Q₁(x,y,z)
# m₂(x,y,z) = M₁₂₁(x,y,z)+M₂₂₂(x,y,z) - Q₂(x,y,z)


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
    prescribe!(elements_w, :E=>E, :ν=>ν, :h=>h, :q=>q)
    prescribe!(elements_φ, :E=>E, :ν=>ν, :h=>h)
    # prescribe!(elements_φ, :E=>E, :ν=>ν, :h=>h, :m₁=>m₁, :m₂=>m₂)
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
    𝑓ʷ = ∫wqdΩ=>elements_w
    𝑓ᵠ = ∫φmdΩ=>elements_φ
    @timeit to "assemble" 𝑎ˢˢ(kˢˢ)
    @timeit to "assemble" 𝑎ˢʷ(kˢʷ)
    @timeit to "assemble" 𝑎ˢᵠ(kˢᵠ)
    @timeit to "assemble" 𝑓ʷ(fʷ)
    # @timeit to "assemble" 𝑓ᵠ(fᵠ)
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
    𝑎ʷʷ = ∫αwwdΓ=>elements_1_w∪elements_2_w∪elements_3_w∪elements_4_w
    𝑎ᵠᵠ = ∫αφφdΓ=>elements_1_φ∪elements_2_φ∪elements_3_φ∪elements_4_φ
    @timeit to "assemble" 𝑎ˢʷ(kˢʷ,fˢ)
    @timeit to "assemble" 𝑎ᵐᵠ(kᵐᵠ,fᵐ)
    @timeit to "assemble" 𝑎ʷʷ(kʷʷ, fʷ)
    @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ, fᵠ)
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
k̃ᵠᵠ = - kᵐᵠ'*(kᵐᵐ\kᵐᵠ)
f̃ᵠ = - kᵐᵠ'*(kᵐᵐ\fᵐ)
@timeit to "solve" d = [k̃ᵠᵠ kᵠʷ kˢᵠ';kᵠʷ' kʷʷ kˢʷ';kˢᵠ kˢʷ kˢˢ]\[fᵠ+f̃ᵠ;fʷ;fˢ]
push!(nodes,:d=>d[2*nᵠ+1:2*nᵠ+nʷ])
# ──────────────────────────────────────────────────────────
@timeit to "calculate error" begin
    @timeit to "get elements" elements_w = getElements(nodes, entities["Ω"], 10)
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"], 10)
    elements_φ, elements_Γ, nodes_φ = Tri3toTRTri3(elements,elements_Γ)
    push!(nodes_φ,:d₁=>d[1:2:2*nᵠ], :d₂=>d[2:2:2*nᵠ], :q₁=>d[2*nᵠ+nʷ+1:2:end], :q₂=>d[2*nᵠ+nʷ+2:2:end])
    prescribe!(elements_w, :E=>E, :ν=>ν, :h=>h, :w=>w, :∂w∂x=>w₁, :∂w∂y=>w₂)
    prescribe!(elements_φ, :E=>E, :ν=>ν, :h=>h, :φ₁=>φ₁, :φ₂=>φ₂, :∂φ₁∂x=>φ₁₁, :∂φ₂∂x=>φ₂₁, :∂φ₁∂y=>φ₁₂, :∂φ₂∂y=>φ₂₂, :Q₁=>Q₁, :Q₂=>Q₂)
    @timeit to "calculate shape functions" set∇𝝭!(elements_w)
    @timeit to "calculate shape functions" set∇𝝭!(elements_φ)

    L₂_w = L₂w(elements_w)
    L₂_φ = L₂φ(elements_φ)
    L₂_Q = L₂Q(elements_φ)
end

println(to)

println("h=$h,Dᵇ=$Dᵇ,Dˢ=$Dˢ,αʷ=$αʷ,αᵠ=$αᵠ,nᵠ=$nᵠ,nʷ=$nʷ,nˢ=$nˢ,nᵐ=$nᵐ")
print("nˢ≤ᵠ:         ")
n_diff = nᵠ-nˢ
n_diff≥0.0 ? println("✓:$n_diff") : println("×:$n_diff")

print("2nˢ≤nʷ+2nᵠ-min(nʷ,n):")
n = floor(0.5*((1+8*nᵠ)^0.5-3))
# n = ceil(0.5*((1+8*nᵠ)^0.5-3))
n = 0.5*(n+2)*(n+3)
n_diff = 0.5*(nʷ+2nᵠ-min(nʷ,n))-nˢ
n_diff≥0.0 ? println("✓:$n_diff") : println("×:$n_diff")
println("nʷ = $nʷ, n = $n")

println("L₂ error of w: ", log10(L₂_w))
println("L₂ error of φ: ", log10(L₂_φ))
println("L₂ error of Q: ", log10(L₂_Q))

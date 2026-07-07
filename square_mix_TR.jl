using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QφdΩ, ∫Q∇wdΩ, ∫QQdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫wVdΓ, ∫QwdΓ, ∫αwwdΓ, ∫αφφdΓ, L₂w, L₂φ, L₂Q, H₁, Hₑ
import ApproxOperator.Heat: ∫uds

using TimerOutputs, WriteVTK
import Gmsh: gmsh

E = 10.92e6
ν = 0.3
h = 1e-5
Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))

αʷ = 1e1*h^2
αᵠ = 1e8
ndiv = 4

# w(x,y,z) = - (x^4*y - x*y^4 + x^3*y^2 - x^2*y^3)/6/Dᵇ + (x^3 + 3*x^2*y - 3*x*y^2 - y^3)/3/Dˢ
# w₁(x,y,z) = - (4*x^3*y - y^4 + 3*x^2*y^2 - 2*x*y^3)/6/Dᵇ + (x^2 + 2*x*y - y^2)/Dˢ
# w₂(x,y,z) = - (x^4 - 4*x*y^3 + 2*x^3*y - 3*x^2*y^2)/6/Dᵇ + (x^2 - 2*x*y - y^2)/Dˢ
# w₁₁(x,y,z) = - (6*x^2*y + 3*x*y^2 - y^3)/3/Dᵇ + 2*(x + y)/Dˢ
# w₂₂(x,y,z) = - (- 6*x*y^2 + x^3 - 3*x^2*y)/3/Dᵇ + 2*(- x - y)/Dˢ
# φ₁(x,y,z) = - (4*x^3*y - y^4 + 3*x^2*y^2 - 2*x*y^3)/6/Dᵇ
# φ₂(x,y,z) = - (x^4 - 4*x*y^3 + 2*x^3*y - 3*x^2*y^2)/6/Dᵇ
# φ₁₁(x,y,z) = - (6*x^2*y + 3*x*y^2 - y^3)/3/Dᵇ
# φ₁₂(x,y,z) = - (2*x^3 - 2*y^3 + 3*x^2*y - 3*x*y^2)/3/Dᵇ
# φ₂₁(x,y,z) = - (2*x^3 - 2*y^3 + 3*x^2*y - 3*x*y^2)/3/Dᵇ
# φ₂₂(x,y,z) = - (- 6*x*y^2 + x^3 - 3*x^2*y)/3/Dᵇ
# φ₁₁₁(x,y,z) = - (4*x*y + y^2)/Dᵇ
# φ₁₁₂(x,y,z) = - (2*x^2 + 2*x*y - y^2)/Dᵇ
# φ₂₂₁(x,y,z) = - (- 2*y^2 + x^2 - 2*x*y)/Dᵇ
# φ₂₂₂(x,y,z) = - (- 4*x*y - x^2)/Dᵇ
# φ₁₂₁(x,y,z) = - (2*x^2 + 2*x*y - y^2)/Dᵇ
# φ₁₂₂(x,y,z) = - (- 2*y^2 + x^2 - 2*x*y)/Dᵇ

# Q₁(x,y,z) = (x^2 + 2*x*y - y^2)
# Q₂(x,y,z) = (x^2 - 2*x*y - y^2)
# q(x,y,z) = 0.0

w(x,y,z) = 1/3*x^3*(x-1)^3*y^3*(y-1)^3-2*h^2/(5*(1-ν))*(y^3*(y-1)^3*x*(x-1)*(5*x^2-5*x+1)+x^3*(x-1)^3*y*(y-1)*(5*y^2-5*y+1))
w₁(x,y,z) = (x-1)^2*x^2*(2*x-1)*(y-1)^3*y^3-2*h^2/(5*(1-ν))*((20*x^3-30*x^2+12*x-1)*(y-1)^3*y^3+3*(x-1)^2*x^2*(2*x-1)*(y-1)*y*(5*y^2-5*y+1))
w₂(x,y,z) = (x-1)^3*x^3*(y-1)^2*y^2*(2*y-1)-2*h^2/(5*(1-ν))*(3*(x-1)*x*(5*x^2-5*x+1)*(y-1)^2*y^2*(2*y-1)+x^3*(x-1)^3*(20*y^3-30*y^2+12*y-1))
φ₁(x,y,z) = y^3*(y-1)^3*x^2*(x-1)^2*(2*x-1)
φ₁₁(x,y,z) = y^3*(y-1)^3*2*(x*(x-1)^2*(2*x-1)+x^2*(x-1)*(2*x-1)+x^2*(x-1)^2)
φ₁₂(x,y,z) = 3*(y^2*(y-1)^3+y^3*(y-1)^2)*x^2*(x-1)^2*(2*x-1)
φ₂(x,y,z) = x^3*(x-1)^3*y^2*(y-1)^2*(2*y-1)
φ₂₁(x,y,z) = 3*(x^2*(x-1)^3+x^3*(x-1)^2)*y^2*(y-1)^2*(2*y-1)
φ₂₂(x,y,z) = x^3*(x-1)^3*2*(y*(y-1)^2*(2*y-1)+y^2*(y-1)*(2*y-1)+y^2*(y-1)^2)
q(x,y,z) = E*h^3/(12*(1-ν^2))*(12*y*(y-1)*(5*x^2-5*x+1)*(2*y^2*(y-1)^2+x*(x-1)*(5*y^2-5*y+1))+12*x*(x-1)*(5*y^2-5*y+1)*(2*x^2*(x-1)^2+y*(y-1)*(5*x^2-5*x+1)))

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

Q₁(x,y,z) = Dˢ*(w₁(x,y,z)-φ₁(x,y,z))
Q₂(x,y,z) = Dˢ*(w₂(x,y,z)-φ₂(x,y,z))
# Q₁₁(x,y,z) = Dˢ*(w₁₁(x,y,z)-φ₁₁(x,y,z))
# Q₂₂(x,y,z) = Dˢ*(w₂₂(x,y,z)-φ₂₂(x,y,z))
# q(x,y,z)=-Q₁₁(x,y,z)-Q₂₂(x,y,z)
# m₁(x,y,z) = M₁₁₁(x,y,z)+M₁₂₂(x,y,z) - Q₁(x,y,z)
# m₂(x,y,z) = M₁₂₁(x,y,z)+M₂₂₂(x,y,z) - Q₂(x,y,z)

integrationOrder = 2

const to = TimerOutput()

gmsh.initialize()
@timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv.msh")
# @timeit to "open msh file" gmsh.open("msh/patchtest_tri3_irregular_1034.msh")
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
    prescribe!(elements_w, :E=>E, :ν=>ν, :h=>h, :q=>q)
    # prescribe!(elements_φ, :E=>E, :ν=>ν, :h=>h, :m₁=>m₁, :m₂=>m₂)
    prescribe!(elements_φ, :E=>E, :ν=>ν, :h=>h)
    @timeit to "calculate shape functions" set∇𝝭!(elements_w)
    @timeit to "calculate shape functions" set∇𝝭!(elements_φ)
    𝑎ᵠᵠ = ∫κκdΩ=>elements_φ
    𝑎ˢᵠ = ∫QφdΩ=>elements_φ
    𝑎ˢˢ = ∫QQdΩ=>elements_φ
    𝑎ˢʷ = ∫Q∇wdΩ=>(elements_φ,elements_w)
    𝑓ʷ = ∫wqdΩ=>elements_w
    # 𝑓ᵠ = ∫φmdΩ=>elements_φ
    @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)
    @timeit to "assemble" 𝑎ˢᵠ(kˢᵠ)
    @timeit to "assemble" 𝑎ˢˢ(kˢˢ)
    @timeit to "assemble" 𝑎ˢʷ(kˢʷ)
    @timeit to "assemble" 𝑓ʷ(fʷ)
    # @timeit to "assemble" 𝑓ᵠ(fᵠ)
end

@timeit to "calculate ∫αwwdΓ" begin
    @timeit to "get elements" elements_1_w = getElements(nodes, entities["Γ¹"], integrationOrder)
    @timeit to "get elements" elements_2_w = getElements(nodes, entities["Γ²"], integrationOrder)
    @timeit to "get elements" elements_3_w = getElements(nodes, entities["Γ³"], integrationOrder)
    @timeit to "get elements" elements_4_w = getElements(nodes, entities["Γ⁴"], integrationOrder)
    elements_1_φ = Seg2toTRTri3(elements_1_w,elements_φ)
    elements_2_φ = Seg2toTRTri3(elements_2_w,elements_φ)
    elements_3_φ = Seg2toTRTri3(elements_3_w,elements_φ)
    elements_4_φ = Seg2toTRTri3(elements_4_w,elements_φ)
    @timeit to "calculate shape functions" set𝝭!(elements_1_w)
    @timeit to "calculate shape functions" set𝝭!(elements_2_w)
    @timeit to "calculate shape functions" set𝝭!(elements_3_w)
    @timeit to "calculate shape functions" set𝝭!(elements_4_w)
    @timeit to "calculate shape functions" set𝝭!(elements_1_φ)
    @timeit to "calculate shape functions" set𝝭!(elements_2_φ)
    @timeit to "calculate shape functions" set𝝭!(elements_3_φ)
    @timeit to "calculate shape functions" set𝝭!(elements_4_φ)
    prescribe!(elements_1_w, :α=>αʷ*Dˢ, :g=>w)
    prescribe!(elements_2_w, :α=>αʷ*Dˢ, :g=>w)
    prescribe!(elements_3_w, :α=>αʷ*Dˢ, :g=>w)
    prescribe!(elements_4_w, :α=>αʷ*Dˢ, :g=>w)
    prescribe!(elements_1_φ, :α=>αᵠ*E, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_2_φ, :α=>αᵠ*E, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_3_φ, :α=>αᵠ*E, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_4_φ, :α=>αᵠ*E, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    𝑎 = ∫QwdΓ => (elements_1_φ∪elements_2_φ∪elements_3_φ∪elements_4_φ,elements_1_w∪elements_2_w∪elements_3_w∪elements_4_w)
    𝑎ʷ = ∫αwwdΓ=>elements_1_w∪elements_2_w∪elements_3_w∪elements_4_w
    𝑎ᵠ = ∫αφφdΓ=>elements_1_φ∪elements_2_φ∪elements_3_φ∪elements_4_φ
    @timeit to "assemble" 𝑎(kˢʷ,fˢ)
    @timeit to "assemble" 𝑎ʷ(kʷʷ,fʷ)
    @timeit to "assemble" 𝑎ᵠ(kᵠᵠ,fᵠ)
end

@timeit to "solve" d = [kᵠᵠ kᵠʷ kˢᵠ';kᵠʷ' kʷʷ kˢʷ';kˢᵠ kˢʷ kˢˢ]\[fᵠ;fʷ;fˢ]
push!(nodes,:d=>d[2*nᵠ+1:2*nᵠ+nʷ])

@timeit to "calculate error" begin
    # prescribe!(elements_Γ, :u=>(x,y,z)->y)
    # u = ∫uds(elements_Γ)
    @timeit to "get elements" elements_w = getElements(nodes, entities["Ω"], 10)
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"], 10)
    elements_φ, elements_Γ, nodes = Tri3toTRTri3(elements,elements_Γ)
    push!(nodes,:d₁=>d[1:2:2*nᵠ], :d₂=>d[2:2:2*nᵠ], :q₁=>d[2*nᵠ+nʷ+1:2:end], :q₂=>d[2*nᵠ+nʷ+2:2:end])
    prescribe!(elements_w, :E=>E, :ν=>ν, :h=>h, :w=>w, :∂w∂x=>w₁, :∂w∂y=>w₂)
    prescribe!(elements_φ, :E=>E, :ν=>ν, :h=>h, :φ₁=>φ₁, :φ₂=>φ₂, :∂φ₁∂x=>φ₁₁, :∂φ₂∂x=>φ₂₁, :∂φ₁∂y=>φ₁₂, :∂φ₂∂y=>φ₂₂, :Q₁=>Q₁, :Q₂=>Q₂)
    @timeit to "calculate shape functions" set∇𝝭!(elements_w)
    @timeit to "calculate shape functions" set∇𝝭!(elements_φ)
    L₂_w = L₂w(elements_w)
    L₂_φ = L₂φ(elements_φ)
    L₂_Q = L₂Q(elements_φ)
    # H₁_u, L₂_u = H₁(elements_w,elements_φ)
    # Hₑ_u = Hₑ(elements_w,elements_φ,elements_φ)
end

gmsh.finalize()

println("h=$h,Dᵇ=$Dᵇ,Dˢ=$Dˢ,nᵠ=$nᵠ,nʷ=$nʷ,nˢ=$nˢ")

print("nʷ≤⌊[nˢ]⌋-1:         ")
n = floor(0.5*((1+8*nˢ)^0.5-3))
n_diff = 0.5*n*(n+1)-nʷ
n_diff≥0.0 ? println("✓:$n_diff") : println("×:$n_diff")

print("2nˢ≤nʷ+2nᵠ-min(nʷ,n):")
n = floor(0.5*((1+8*nᵠ)^0.5-3))
# n = ceil(0.5*((1+8*nᵠ)^0.5-3))
n = 0.5*(n+2)*(n+3)
n_diff = 0.5*(nʷ+2nᵠ-min(nʷ,n))-nˢ
n_diff≥0.0 ? println("✓:$n_diff") : println("×:$n_diff")
println("nʷ = $nʷ, n = $n")

# h = log10(1/ndiv)
L₂_w = log10(L₂_w)
L₂_φ = log10(L₂_φ)
L₂_Q = log10(L₂_Q)
# L₂_u = log10(L₂_u)
# H₁_u = log10(H₁_u)
# Hₑ_u = log10(Hₑ_u)

# println("h=$h, L₂_w=$L₂_w, L₂_φ=$L₂_φ, L₂_Q=$L₂_Q")
println("$L₂_w, $L₂_φ, $L₂_Q")
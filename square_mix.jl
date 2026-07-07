using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, ∫wVdΓ, ∫φMdΓ, L₂w, L₂φ, L₂Q

using TimerOutputs, LinearAlgebra
import Gmsh: gmsh

E = 10.92e6
ν = 0.3
h = 1e-0
Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))

αʷ = 1e1*h^2
ndiv = 4
αᵠ = 1e8

w(x,y,z) = - (x^4*y - x*y^4 + x^3*y^2 - x^2*y^3)/6/Dᵇ + (x^3 + 3*x^2*y - 3*x*y^2 - y^3)/3/Dˢ
w₁(x,y,z) = - (4*x^3*y - y^4 + 3*x^2*y^2 - 2*x*y^3)/6/Dᵇ + (x^2 + 2*x*y - y^2)/Dˢ
w₂(x,y,z) = - (x^4 - 4*x*y^3 + 2*x^3*y - 3*x^2*y^2)/6/Dᵇ + (x^2 - 2*x*y - y^2)/Dˢ
w₁₁(x,y,z) = - (6*x^2*y + 3*x*y^2 - y^3)/3/Dᵇ + 2*(x + y)/Dˢ
w₂₂(x,y,z) = - (- 6*x*y^2 + x^3 - 3*x^2*y)/3/Dᵇ + 2*(- x - y)/Dˢ
φ₁(x,y,z) = - (4*x^3*y - y^4 + 3*x^2*y^2 - 2*x*y^3)/6/Dᵇ
φ₂(x,y,z) = - (x^4 - 4*x*y^3 + 2*x^3*y - 3*x^2*y^2)/6/Dᵇ
φ₁₁(x,y,z) = - (6*x^2*y + 3*x*y^2 - y^3)/3/Dᵇ
φ₁₂(x,y,z) = - (2*x^3 - 2*y^3 + 3*x^2*y - 3*x*y^2)/3/Dᵇ
φ₂₁(x,y,z) = - (2*x^3 - 2*y^3 + 3*x^2*y - 3*x*y^2)/3/Dᵇ
φ₂₂(x,y,z) = - (- 6*x*y^2 + x^3 - 3*x^2*y)/3/Dᵇ
φ₁₁₁(x,y,z) = - (4*x*y + y^2)/Dᵇ
φ₁₁₂(x,y,z) = - (2*x^2 + 2*x*y - y^2)/Dᵇ
φ₂₂₁(x,y,z) = - (- 2*y^2 + x^2 - 2*x*y)/Dᵇ
φ₂₂₂(x,y,z) = - (- 4*x*y - x^2)/Dᵇ
φ₁₂₁(x,y,z) = - (2*x^2 + 2*x*y - y^2)/Dᵇ
φ₁₂₂(x,y,z) = - (- 2*y^2 + x^2 - 2*x*y)/Dᵇ

Q₁(x,y,z) = (x^2 + 2*x*y - y^2)
Q₂(x,y,z) = (x^2 - 2*x*y - y^2)
q(x,y,z) = 0.0

# w(x,y,z) = 1/3*x^3*(x-1)^3*y^3*(y-1)^3-2*h^2/(5*(1-ν))*(y^3*(y-1)^3*x*(x-1)*(5*x^2-5*x+1)+x^3*(x-1)^3*y*(y-1)*(5*y^2-5*y+1))
# w₁(x,y,z) = (x-1)^2*x^2*(2*x-1)*(y-1)^3*y^3-2*h^2/(5*(1-ν))*((20*x^3-30*x^2+12*x-1)*(y-1)^3*y^3+3*(x-1)^2*x^2*(2*x-1)*(y-1)*y*(5*y^2-5*y+1))
# w₂(x,y,z) = (x-1)^3*x^3*(y-1)^2*y^2*(2*y-1)-2*h^2/(5*(1-ν))*(3*(x-1)*x*(5*x^2-5*x+1)*(y-1)^2*y^2*(2*y-1)+x^3*(x-1)^3*(20*y^3-30*y^2+12*y-1))
# φ₁(x,y,z) = y^3*(y-1)^3*x^2*(x-1)^2*(2*x-1)
# φ₂(x,y,z) = x^3*(x-1)^3*y^2*(y-1)^2*(2*y-1)
# q(x,y,z) = E*h^3/(12*(1-ν^2))*(12*y*(y-1)*(5*x^2-5*x+1)*(2*y^2*(y-1)^2+x*(x-1)*(5*y^2-5*y+1))+12*x*(x-1)*(5*y^2-5*y+1)*(2*x^2*(x-1)^2+y*(y-1)*(5*x^2-5*x+1)))

# Q₁(x,y,z) = Dˢ*(w₁(x,y,z)-φ₁(x,y,z)) - (x^2 + 2*x*y - y^2)
# Q₂(x,y,z) = Dˢ*(w₂(x,y,z)-φ₂(x,y,z)) - (x^2 - 2*x*y - y^2)

# w(x,y,z) = 1.0+x+y
# w₁(x,y,z) = 1.0
# w₂(x,y,z) = 1.0

# φ₁(x,y,z) = 1.0+x+y
# φ₂(x,y,z) = 1.0+x+y

# Q₁(x,y,z) = Dˢ*(w₁(x,y,z)-φ₁(x,y,z))
# Q₂(x,y,z) = Dˢ*(w₂(x,y,z)-φ₂(x,y,z))

# q(x,y,z) = 0.0
const to = TimerOutput()

gmsh.initialize()
# @timeit to "open msh file" gmsh.open("msh/patchtest_tri3_irregular_1034.msh")
@timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv.msh")
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()

nʷ = length(nodes)
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
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"], integrationOrder)
    @timeit to "get elements" elements_Γ = getElements(nodes, entities["Γ"], integrationOrder, normal=true)
    prescribe!(elements, :E=>E, :ν=>ν, :h=>h, :q=>q)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    @timeit to "calculate shape functions" set𝝭!(elements_Γ)
    𝑎ᵠᵠ = ∫κκdΩ=>elements
    𝑎ˢᵠ = ∫QφdΩ=>elements
    𝑎ˢˢ = ∫QQdΩ=>elements
    𝑎ˢʷ = [
        ∫∇QwdΩ=>elements,
        ∫QwdΓ=>elements_Γ,
    ]
    𝑓ʷ = ∫wqdΩ=>elements
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
    prescribe!(elements_1, :α=>αᵠ*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_2, :α=>αᵠ*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_3, :α=>αᵠ*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_4, :α=>αᵠ*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    @timeit to "calculate shape functions" set𝝭!(elements_4)
    𝑎ᵠ = ∫αφφdΓ=>elements_1∪elements_2∪elements_3∪elements_4
    @timeit to "assemble" 𝑎ᵠ(kᵠᵠ,fᵠ)
    𝑎ˢ = ∫QwdΓ=>elements_1∪elements_2∪elements_3∪elements_4
    @timeit to "assemble" 𝑎ˢ(kˢʷ,fˢ)
    prescribe!(elements_1, :α=>αʷ*Dˢ, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_2, :α=>αʷ*Dˢ, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_3, :α=>αʷ*Dˢ, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_4, :α=>αʷ*Dˢ, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    𝑎ʷ = ∫αwwdΓ=>elements_1∪elements_2∪elements_3∪elements_4
    @timeit to "assemble" 𝑎ʷ(kʷʷ,fʷ)
end

dᵠ = zeros(2*nᵠ)
dˢ = zeros(2*nˢ)
dʷ = zeros(nʷ)
for node in nodes
    x = node.x
    y = node.y
    z = node.z
    dᵠ[2*node.𝐼-1] = φ₁(x,y,z)
    dᵠ[2*node.𝐼]   = φ₂(x,y,z)
    dˢ[2*node.𝐼-1] = Q₁(x,y,z)
    dˢ[2*node.𝐼]   = Q₂(x,y,z)
    dʷ[node.𝐼] = w(x,y,z)
end
# println(kˢˢ*dˢ)
# println(kˢʷ*dʷ)
# println(kˢˢ*dˢ + kˢʷ*dʷ)
# println(kˢʷ*ones(nʷ).-fˢ)
# println(kᵠᵠ*dᵠ + kˢᵠ'*dˢ - fᵠ)
# println(kˢˢ*dˢ + kˢᵠ*dᵠ + kˢʷ*dʷ - fˢ)
# println(kˢʷ'*dˢ + kʷʷ*dʷ - fʷ)
# println(norm(kᵠᵠ*dᵠ + kˢᵠ'*dˢ - fᵠ))
# println(norm(kˢˢ*dˢ + kˢᵠ*dᵠ + kˢʷ*dʷ - fˢ))
# println(norm(kˢʷ'*dˢ + kʷʷ*dʷ - fʷ))
# println(kᵠᵠ*dᵠ + kˢᵠ'*dˢ - fᵠ)
# println(kᵠᵠ*dᵠ-fᵠ)
# println(kᵠᵠ*dᵠ)
# println(fᵠ)
# println(kˢʷ'*dˢ)
# println(kˢᵠ*dᵠ)
# println(kˢʷ*dʷ)
# println(kˢˢ*dˢ)
# println(kˢˢ*dˢ + kˢʷ*dʷ)

@timeit to "solve" d = [kᵠᵠ kᵠʷ kˢᵠ';kᵠʷ' kʷʷ kˢʷ';kˢᵠ kˢʷ kˢˢ]\[fᵠ;fʷ;fˢ]
push!(nodes,:d=>d[2*nᵠ+1:2*nᵠ+nʷ], :d₁=>d[1:2:2*nᵠ], :d₂=>d[2:2:2*nᵠ], :q₁=>d[2*nᵠ+nʷ+1:2:end], :q₂=>d[2*nᵠ+nʷ+2:2:end])

@timeit to "calculate error" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"], 10)
    prescribe!(elements, :E=>E, :ν=>ν, :h=>h, :w=>w, :φ₁=>φ₁, :φ₂=>φ₂, :Q₁=>Q₁, :Q₂=>Q₂)
    @timeit to "calculate shape functions" set𝝭!(elements)
    L₂_w = L₂w(elements)
    L₂_φ = L₂φ(elements)
    L₂_Q = L₂Q(elements)
end

gmsh.finalize()

# points = zeros(3, nʷ)
# for node in nodes
#     I = node.𝐼
#     points[1,I] = node.x
#     points[2,I] = node.y
#     points[3,I] = node.z
# end
# # cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, [node.𝐼 for node in elm.𝓒]) for elm in elements]
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [node.𝐼 for node in elm.𝓒]) for elm in elements]
# vtk_grid("vtk/square.vtu", points, cells) do vtk
#     vtk["Q₁"] = [node.q₁ for node in nodes]
#     vtk["Q₂"] = [node.q₂ for node in nodes]
#     vtk["Q̄₁"] = [Q₁(node.x,node.y,node.z) for node in nodes]
#     vtk["Q̄₂"] = [Q₂(node.x,node.y,node.z) for node in nodes]
# end

# println(to)

# println("L₂ error of w: ", L₂_w)
# println("L₂ error of φ: ", L₂_φ)
# println("L₂ error of Q: ", L₂_Q)

logL₂w = log10(L₂_w)
logL₂φ = log10(L₂_φ)
logL₂Q = log10(L₂_Q)
println("$logL₂w, $logL₂φ, $logL₂Q")


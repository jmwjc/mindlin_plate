using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QQdΩ, ∫Q∇wdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, ∫wVdΓ, ∫φMdΓ, L₂w, L₂φ, L₂Q

using TimerOutputs, WriteVTK, XLSX 
import Gmsh: gmsh

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

Q₁(x,y,z) = Dˢ*(w₁(x,y,z)-φ₁(x,y,z))
Q₂(x,y,z) = Dˢ*(w₂(x,y,z)-φ₂(x,y,z))

const to = TimerOutput()

gmsh.initialize()

integrationOrder = 2
type_w = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_φ = :tri3
type_q = :tri3
ndiv_w = 4

@timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv_w.msh")
@timeit to "get nodes" nodes_w = get𝑿ᵢ()
xʷ = nodes_w.x
yʷ = nodes_w.y
zʷ = nodes_w.z
nʷ = length(nodes_w)
sp = RegularGrid(xʷ,yʷ,zʷ,n = 3,γ = 5)
s = 1/ndiv_w
s₁ = 1.5*s*ones(nʷ)
s₂ = 1.5*s*ones(nʷ)
s₃ = 1.5*s*ones(nʷ)
push!(nodes_w,:s₁=>s₁,:s₂=>s₂,:s₃=>s₃)

@timeit to "open msh file" gmsh.open("msh/patchtest_$type_φ"*"_$ndiv.msh")
@timeit to "get nodes" nodes_φ = get𝑿ᵢ()
@timeit to "get entities" entities = getPhysicalGroups()

@timeit to "calculate main elements" begin
    @timeit to "get elements" elements_φ = getElements(nodes_φ, entities["Ω"], integrationOrder)
    @timeit to "get elements" elements_w = getElements(nodes_w, entities["Ω"], eval(type_w), integrationOrder, sp)
    @timeit to "get elements" elements_q = getPiecewiseElements(entities["Ω"], eval(type_q), integrationOrder)
end
nₑ = length(elements_φ)
nᵠ = length(nodes_φ)
nˢ = nₑ*ApproxOperator.get𝑛𝑝(elements_q[1])
kʷʷ = zeros(nʷ,nʷ)
kˢʷ = zeros(2*nˢ,nʷ)
kᵠʷ = zeros(2*nᵠ,nʷ)
fʷ = zeros(nʷ)
kᵠᵠ = zeros(2*nᵠ,2*nᵠ)
kˢˢ = zeros(2*nˢ,2*nˢ)
kˢᵠ = zeros(2*nˢ,2*nᵠ)
fᵠ = zeros(2*nᵠ)
fˢ = zeros(2*nˢ)

@timeit to "calculate ∫κκdΩ" begin
    @timeit to "get elements" elements_w_Γ = getElements(nodes_w, entities["Γ"], eval(type_w), integrationOrder, sp, normal=true)
    @timeit to "get elements" elements_q_Γ = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], eval(type_q), integrationOrder)
    @timeit to "get elements" elements_φ_Γ = getElements(nodes_φ, entities["Γ"], integrationOrder, normal=true)
    prescribe!(elements_φ, :E=>E, :ν=>ν, :h=>h)
    prescribe!(elements_q, :E=>E, :ν=>ν, :h=>h)
    prescribe!(elements_w, :E=>E, :ν=>ν, :h=>h, :q=>q)
    @timeit to "calculate shape functions" set∇𝝭!(elements_φ)
    @timeit to "calculate shape functions" set∇𝝭!(elements_q)
    @timeit to "calculate shape functions" set∇𝝭!(elements_w)
    @timeit to "calculate shape functions" set𝝭!(elements_w_Γ)
    @timeit to "calculate shape functions" set𝝭!(elements_q_Γ)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_Γ)
    𝑎ᵠᵠ = ∫κκdΩ=>elements_φ
    𝑎ˢᵠ = ∫QφdΩ=>(elements_q,elements_φ)
    𝑎ˢˢ = ∫QQdΩ=>elements_q
    # 𝑎ˢʷ = ∫Q∇wdΩ=>(elements_q,elements_w)
    𝑎ˢʷ = [
        ∫∇QwdΩ=>(elements_q,elements_w),
        ∫QwdΓ=>(elements_q_Γ,elements_w_Γ),
    ]
    # 𝑓ᵠ = ∫φmdΩ=>elements_φ
    𝑓ʷ = ∫wqdΩ=>elements_w
    @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)
    @timeit to "assemble" 𝑎ˢˢ(kˢˢ)
    @timeit to "assemble" 𝑎ˢᵠ(kˢᵠ)
    @timeit to "assemble" 𝑎ˢʷ(kˢʷ)
    @timeit to "assemble" 𝑓ʷ(fʷ)
    # @timeit to "assemble" 𝑓ᵠ(fᵠ)
end

@timeit to "calculate ∫αφφdΓ" begin
    @timeit to "get elements" elements_1 = getElements(nodes_φ, entities["Γ¹"], integrationOrder, normal=true)
    @timeit to "get elements" elements_2 = getElements(nodes_φ, entities["Γ²"], integrationOrder, normal=true)
    @timeit to "get elements" elements_3 = getElements(nodes_φ, entities["Γ³"], integrationOrder, normal=true)
    @timeit to "get elements" elements_4 = getElements(nodes_φ, entities["Γ⁴"], integrationOrder, normal=true)
    prescribe!(elements_1, :α=>1e8*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_2, :α=>1e8*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_3, :α=>1e8*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_4, :α=>1e8*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    @timeit to "calculate shape functions" set𝝭!(elements_4)
    𝑎ᵠ = ∫αφφdΓ=>elements_1∪elements_2∪elements_3∪elements_4
    @timeit to "assemble" 𝑎ᵠ(kᵠᵠ,fᵠ)
end

@timeit to "calculate ∫αwwdΓ ∫QwdΓ" begin
    @timeit to "get elements" elements_q_1 = getElements(entities["Γ¹"], entities["Γ"], elements_q_Γ)
    @timeit to "get elements" elements_q_2 = getElements(entities["Γ²"], entities["Γ"], elements_q_Γ)
    @timeit to "get elements" elements_q_3 = getElements(entities["Γ³"], entities["Γ"], elements_q_Γ)
    @timeit to "get elements" elements_q_4 = getElements(entities["Γ⁴"], entities["Γ"], elements_q_Γ)
    @timeit to "get elements" elements_w_1 = getElements(nodes_w, entities["Γ¹"], eval(type_w), integrationOrder, sp, normal=true)
    @timeit to "get elements" elements_w_2 = getElements(nodes_w, entities["Γ²"], eval(type_w), integrationOrder, sp, normal=true)
    @timeit to "get elements" elements_w_3 = getElements(nodes_w, entities["Γ³"], eval(type_w), integrationOrder, sp, normal=true)
    @timeit to "get elements" elements_w_4 = getElements(nodes_w, entities["Γ⁴"], eval(type_w), integrationOrder, sp, normal=true)
    prescribe!(elements_w_1, :α=>1e8*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_w_2, :α=>1e8*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_w_3, :α=>1e8*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_w_4, :α=>1e8*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    @timeit to "calculate shape functions" set𝝭!(elements_w_1)
    @timeit to "calculate shape functions" set𝝭!(elements_w_2)
    @timeit to "calculate shape functions" set𝝭!(elements_w_3)
    @timeit to "calculate shape functions" set𝝭!(elements_w_4)
    𝑎ˢ = ∫QwdΓ=>(elements_q_1∪elements_q_2∪elements_q_3∪elements_q_4,elements_w_1∪elements_w_2∪elements_w_3∪elements_w_4)
    @timeit to "assemble" 𝑎ˢ(kˢʷ,fˢ)
end

@timeit to "calculate error" begin
    @timeit to "get elements" elements_φ = getElements(nodes_φ, entities["Ω"], 10)
    @timeit to "get elements" elements_w = getElements(nodes_w, entities["Ω"], eval(type_w), 10, sp)
    @timeit to "get elements" elements_q = getPiecewiseElements(entities["Ω"], eval(type_q), 10)
    prescribe!(elements_φ, :E=>E, :ν=>ν, :h=>h, :φ₁=>φ₁, :φ₂=>φ₂)
    @timeit to "calculate shape functions" set𝝭!(elements_φ)
    prescribe!(elements_w, :E=>E, :ν=>ν, :h=>h, :w=>w)
    @timeit to "calculate shape functions" set𝝭!(elements_w)
    prescribe!(elements_q, :E=>E, :ν=>ν, :h=>h, :Q₁=>Q₁, :Q₂=>Q₂)
    @timeit to "calculate shape functions" set𝝭!(elements_q)
end

@timeit to "solve" d = [kᵠᵠ kᵠʷ kˢᵠ';kᵠʷ' kʷʷ kˢʷ';kˢᵠ kˢʷ kˢˢ]\[fᵠ;fʷ;fˢ]
# println([kᵠᵠ kᵠʷ kˢᵠ';kᵠʷ' kʷʷ kˢʷ';kˢᵠ kˢʷ kˢˢ]*d .- [fᵠ;fʷ;fˢ])
nodes_q = 𝑿ᵢ[]
for elm in elements_q
    for node in elm.𝓒
        push!(nodes_q, node)
    end
end
push!(nodes_φ,:d₁=>d[1:2:2*nᵠ], :d₂=>d[2:2:2*nᵠ])
push!(nodes_w,:d=>d[2*nᵠ+1:2*nᵠ+nʷ])
push!(nodes_q,:q₁=>d[2*nᵠ+nʷ+1:2:end], :q₂=>d[2*nᵠ+nʷ+2:2:end])

@timeit to "calculate error" begin
    L₂_w = L₂w(elements_w)
    L₂_φ = L₂φ(elements_φ)
    L₂_Q = L₂Q(elements_q)
end

# dᵠ = zeros(2*nᵠ)
# dˢ = zeros(2*nˢ)
# dʷ = zeros(nʷ)
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

# for i in 1:nₑ
#     dˢ[6*i-5:2:6*i] = [Q₁(0,0,0),Q₁(1,0,0)-Q₁(0,0,0),Q₁(0,1,0)-Q₁(0,0,0)]
#     dˢ[6*i-4:2:6*i] = [Q₂(0,0,0),Q₂(1,0,0)-Q₂(0,0,0),Q₂(0,1,0)-Q₂(0,0,0)]
# end

# println(kᵠᵠ*dᵠ+kˢᵠ'*dˢ - fᵠ)
# println(norm(kᵠᵠ*dᵠ+kᵠʷ*dʷ+kˢᵠ'*dˢ - fᵠ))
# println(kˢˢ*dˢ)
# println(kˢʷ*dʷ)
# println(kˢˢ*dˢ + kˢʷ*dʷ)
# println(kˢʷ*ones(nʷ).-fˢ)
# println(kᵠᵠ*dᵠ + kˢᵠ'*dˢ - fᵠ)
# println(kˢˢ*dˢ + kˢᵠ*dᵠ + kˢʷ*dʷ - fˢ)
# println(kˢʷ'*dˢ + kʷʷ*dʷ - fʷ)
# println(kˢᵠ*dᵠ)
# println(kˢʷ*dʷ)
# println(kᵠʷ*dʷ)
# err = kˢʷ*dʷ
# println(dˢ'*kˢʷ)
# println(kˢˢ*dˢ)
# println(kˢˢ*dˢ + kˢʷ*dʷ)


# points = zeros(3, nˢ)
# for node in nodes_q
#     I = node.𝐼
#     points[1,I] = node.x
#     points[2,I] = node.y
#     points[3,I] = node.z
# end
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, [node.𝐼 for node in elm.𝓒]) for elm in elements]
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [node.𝐼 for node in elm.𝓒]) for elm in elements_q]
# vtk_grid("vtk/square.vtu", points, cells) do vtk
#     vtk["Q₁"] = [node.q₁ for node in nodes_q]
#     vtk["Q₂"] = [node.q₂ for node in nodes_q]
#     vtk["Q̄₁"] = [Q₁(node.x,node.y,node.z) for node in nodes_q]
#     vtk["Q̄₂"] = [Q₂(node.x,node.y,node.z) for node in nodes_q]
# end

# println(to)

println("L₂ error of w: ", L₂_w)
println("L₂ error of φ: ", L₂_φ)
println("L₂ error of Q: ", L₂_Q)
# ──────────────────────────────────────────────────────────
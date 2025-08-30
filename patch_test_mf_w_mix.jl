using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QQdΩ, ∫Q∇wdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, ∫wVdΓ, ∫φMdΓ, L₂, L₂φ, L₂Q

using TimerOutputs, WriteVTK, XLSX 
import Gmsh: gmsh

E = 1.0
ν = 0.3
h = 1e-8
Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))

r = 1
w(x,y,z) = (x+y)^r
w₁(x,y,z) = r*(x+y)^abs(r-1)
w₂(x,y,z) = r*(x+y)^abs(r-1)
w₁₁(x,y,z) = r*(r-1)*(x+y)^abs(r-2)
w₂₂(x,y,z) = r*(r-1)*(x+y)^abs(r-2)
φ₁(x,y,z) = r*(x+y)^abs(r-1)
φ₂(x,y,z) = r*(x+y)^abs(r-1)
φ₁₁(x,y,z)  = r*(r-1)*(x+y)^abs(r-2)
φ₁₂(x,y,z)  = r*(r-1)*(x+y)^abs(r-2)
φ₂₁(x,y,z)  = r*(r-1)*(x+y)^abs(r-2)
φ₂₂(x,y,z)  = r*(r-1)*(x+y)^abs(r-2)
φ₁₁₁(x,y,z)  = r*(r-1)*(r-2)*(x+y)^abs(r-3)
φ₁₁₂(x,y,z)  = r*(r-1)*(r-2)*(x+y)^abs(r-3)
φ₂₂₁(x,y,z)  = r*(r-1)*(r-2)*(x+y)^abs(r-3)
φ₂₂₂(x,y,z)  = r*(r-1)*(r-2)*(x+y)^abs(r-3)
φ₁₂₁(x,y,z)  = r*(r-1)*(r-2)*(x+y)^abs(r-3)
φ₁₂₂(x,y,z)  = r*(r-1)*(r-2)*(x+y)^abs(r-3)

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

# w(x,y,z) = 2.0
# w₁(x,y,z) = 0.0
# w₂(x,y,z) = 0.0
# w₁₁(x,y,z) = 0.0
# w₂₂(x,y,z) = 0.0
# φ₁(x,y,z) = 1.0
# φ₂(x,y,z) = 1.0
# φ₁₁(x,y,z)  = 0.0
# φ₁₂(x,y,z)  = 0.0
# φ₂₁(x,y,z)  = 0.0
# φ₂₂(x,y,z)  = 0.0
# φ₁₁₁(x,y,z)  = 0.0
# φ₁₁₂(x,y,z)  = 0.0
# φ₂₂₁(x,y,z)  = 0.0
# φ₂₂₂(x,y,z)  = 0.0
# φ₁₂₁(x,y,z)  = 0.0
# φ₁₂₂(x,y,z)  = 0.0

M₁₁(x,y,z)= -Dᵇ*(φ₁₁(x,y,z)+ν*φ₂₂(x,y,z))
M₁₂(x,y,z)= -Dᵇ*(1-ν)*0.5*(φ₁₂(x,y,z)+φ₂₁(x,y,z))
M₂₂(x,y,z)= -Dᵇ*(ν*φ₁₁(x,y,z)+φ₂₂(x,y,z))
M₁₁₁(x,y,z)= -Dᵇ*(φ₁₁₁(x,y,z)+ν*φ₂₂₁(x,y,z))
M₁₂₂(x,y,z)= -Dᵇ*(1-ν)*φ₁₂₂(x,y,z)
M₁₂₁(x,y,z)= -Dᵇ*(1-ν)*φ₁₂₁(x,y,z)
M₂₂₂(x,y,z)= -Dᵇ*(ν*φ₁₁₂(x,y,z)+φ₂₂₂(x,y,z))

Q₁(x,y,z) = Dˢ*(w₁(x,y,z)-φ₁(x,y,z))
Q₂(x,y,z) = Dˢ*(w₂(x,y,z)-φ₂(x,y,z))
Q₁₁(x,y,z) = Dˢ*(w₁₁(x,y,z)-φ₁₁(x,y,z))
Q₂₂(x,y,z) = Dˢ*(w₂₂(x,y,z)-φ₂₂(x,y,z))
q(x,y,z)=-Q₁₁(x,y,z)-Q₂₂(x,y,z)
m₁(x,y,z) = M₁₁₁(x,y,z)+M₁₂₂(x,y,z) - Q₁(x,y,z)
m₂(x,y,z) = M₁₂₁(x,y,z)+M₂₂₂(x,y,z) - Q₂(x,y,z)


const to = TimerOutput()

gmsh.initialize()
# @timeit to "open msh file" gmsh.open("msh/patchtest_3.msh")
# @timeit to "get nodes" nodes_s = get𝑿ᵢ()

integrationOrder = 4
type_w = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_φ = :tri3
# type_q = :(PiecewisePolynomial{:Linear2D})
type_q = :(PiecewisePolynomial{:Quadratic2D})
ndiv = 8
ndiv_w = 4
# XLSX.openxlsx("xls/patchtest.xlsx", mode="w") do xf
# for ndiv_w = 2:42
# row = ndiv_w
# ──────────────────────────────────────────────────────────
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
nᵛ = nₑ*ApproxOperator.get𝑛𝑝(elements_q[1])
kʷʷ = zeros(nʷ,nʷ)
kᵛʷ = zeros(2*nᵛ,nʷ)
kᵠʷ = zeros(2*nᵠ,nʷ)
fʷ = zeros(nʷ)
kᵠᵠ = zeros(2*nᵠ,2*nᵠ)
kᵛᵛ = zeros(2*nᵛ,2*nᵛ)
kᵛᵠ = zeros(2*nᵛ,2*nᵠ)
fᵠ = zeros(2*nᵠ)
fᵛ = zeros(2*nᵛ)

@timeit to "calculate ∫κκdΩ" begin
    @timeit to "get elements" elements_w_Γ = getElements(nodes_w, entities["Γ"], eval(type_w), integrationOrder, sp, normal=true)
    @timeit to "get elements" elements_q_Γ = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], eval(type_q), integrationOrder)
    @timeit to "get elements" elements_φ_Γ = getElements(nodes_φ, entities["Γ"], integrationOrder, normal=true)
    prescribe!(elements_φ, :E=>E, :ν=>ν, :h=>h, :m₁=>m₁, :m₂=>m₂)
    prescribe!(elements_q, :E=>E, :ν=>ν, :h=>h)
    prescribe!(elements_w, :E=>E, :ν=>ν, :h=>h, :m₁=>m₁, :m₂=>m₂, :q=>q)
    @timeit to "calculate shape functions" set∇𝝭!(elements_φ)
    @timeit to "calculate shape functions" set∇𝝭!(elements_q)
    @timeit to "calculate shape functions" set∇𝝭!(elements_w)
    @timeit to "calculate shape functions" set𝝭!(elements_w_Γ)
    @timeit to "calculate shape functions" set𝝭!(elements_q_Γ)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_Γ)
    𝑎ᵠᵠ = ∫κκdΩ=>elements_φ
    𝑎ᵛᵠ = ∫QφdΩ=>(elements_q,elements_φ)
    𝑎ᵛᵛ = ∫QQdΩ=>elements_q
    # 𝑎ᵛʷ = ∫Q∇wdΩ=>(elements_q,elements_w)
    𝑎ᵛʷ = [
        ∫∇QwdΩ=>(elements_q,elements_w),
        ∫QwdΓ=>(elements_q_Γ,elements_w_Γ),
    ]
    𝑓ᵠ = ∫φmdΩ=>elements_φ
    𝑓ʷ = ∫wqdΩ=>elements_w
    @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)
    @timeit to "assemble" 𝑎ᵛᵛ(kᵛᵛ)
    @timeit to "assemble" 𝑎ᵛᵠ(kᵛᵠ)
    @timeit to "assemble" 𝑎ᵛʷ(kᵛʷ)
    @timeit to "assemble" 𝑓ʷ(fʷ)
    @timeit to "assemble" 𝑓ᵠ(fᵠ)
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
    𝑎ᵛ = ∫QwdΓ=>(elements_q_1∪elements_q_2∪elements_q_3∪elements_q_4,elements_w_1∪elements_w_2∪elements_w_3∪elements_w_4)
    @timeit to "assemble" 𝑎ᵛ(kᵛʷ,fᵛ)
end

@timeit to "calculate error" begin
    @timeit to "get elements" elements_φ = getElements(nodes_φ, entities["Ω"], 10)
    @timeit to "get elements" elements_w = getElements(nodes_w, entities["Ω"], eval(type_w), 10, sp)
    @timeit to "get elements" elements_q = getPiecewiseElements(entities["Ω"], eval(type_q), 10)
    prescribe!(elements_φ, :E=>E, :ν=>ν, :h=>h, :φ₁=>φ₁, :φ₂=>φ₂)
    @timeit to "calculate shape functions" set𝝭!(elements_φ)
    prescribe!(elements_w, :E=>E, :ν=>ν, :h=>h, :u=>w)
    @timeit to "calculate shape functions" set𝝭!(elements_w)
    prescribe!(elements_q, :E=>E, :ν=>ν, :h=>h, :Q₁=>Q₁, :Q₂=>Q₂)
    @timeit to "calculate shape functions" set𝝭!(elements_q)
end

@timeit to "solve" d = [kᵠᵠ kᵠʷ kᵛᵠ';kᵠʷ' kʷʷ kᵛʷ';kᵛᵠ kᵛʷ kᵛᵛ]\[fᵠ;fʷ;fᵛ]
# println([kᵠᵠ kᵠʷ kᵛᵠ';kᵠʷ' kʷʷ kᵛʷ';kᵛᵠ kᵛʷ kᵛᵛ]*d .- [fᵠ;fʷ;fᵛ])
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
    L₂_w = L₂(elements_w)
    L₂_φ = L₂φ(elements_φ)
    L₂_Q = L₂Q(elements_q)
end

# dᵠ = zeros(2*nᵠ)
# dᵛ = zeros(2*nᵛ)
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
#     dᵛ[6*i-5:2:6*i] = [Q₁(0,0,0),Q₁(1,0,0)-Q₁(0,0,0),Q₁(0,1,0)-Q₁(0,0,0)]
#     dᵛ[6*i-4:2:6*i] = [Q₂(0,0,0),Q₂(1,0,0)-Q₂(0,0,0),Q₂(0,1,0)-Q₂(0,0,0)]
# end

# println(kᵠᵠ*dᵠ+kᵛᵠ'*dᵛ - fᵠ)
# println(norm(kᵠᵠ*dᵠ+kᵠʷ*dʷ+kᵛᵠ'*dᵛ - fᵠ))
# println(kᵛᵛ*dᵛ)
# println(kᵛʷ*dʷ)
# println(kᵛᵛ*dᵛ + kᵛʷ*dʷ)
# println(kᵛʷ*ones(nʷ).-fᵛ)
# println(kᵠᵠ*dᵠ + kᵛᵠ'*dᵛ - fᵠ)
# println(kᵛᵛ*dᵛ + kᵛᵠ*dᵠ + kᵛʷ*dʷ - fᵛ)
# println(kᵛʷ'*dᵛ + kʷʷ*dʷ - fʷ)
# println(kᵛᵠ*dᵠ)
# println(kᵛʷ*dʷ)
# println(kᵠʷ*dʷ)
# err = kᵛʷ*dʷ
# println(dᵛ'*kᵛʷ)
# println(kᵛᵛ*dᵛ)
# println(kᵛᵛ*dᵛ + kᵛʷ*dʷ)


# points = zeros(3, nᵛ)
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

println(to)

println("L₂ error of w: ", L₂_w)
println("L₂ error of φ: ", L₂_φ)
println("L₂ error of Q: ", L₂_Q)
# ──────────────────────────────────────────────────────────
#     sheet = xf[1]
#     XLSX.rename!(sheet, "new_sheet")
#     sheet["A1"] = "type w"
#     sheet["B1"] = "nʷ"
#     sheet["C1"] = "type φ"
#     sheet["D1"] = "nᵠ"
#     sheet["E1"] = "type Q"
#     sheet["F1"] = "nᵛ"
#     sheet["G1"] = "L₂w"
#     sheet["H1"] = "L₂φ"
#     sheet["I1"] = "L₂Q"
#     sheet["A$row"] = "$type_w"
#     sheet["B$row"] = nʷ
#     sheet["C$row"] = "$type_φ"
#     sheet["D$row"] = nᵠ
#     sheet["E$row"] = "$type_q"
#     sheet["F$row"] = nᵛ
#     sheet["G$row"] = log10(L₂_w)
#     sheet["H$row"] = log10(L₂_φ)
#     sheet["I$row"] = log10(L₂_Q)
# end
# end
gmsh.finalize()



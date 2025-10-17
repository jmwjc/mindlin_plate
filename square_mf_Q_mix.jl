using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QQdΩ, ∫Q∇wdΩ, ∫QwdΓ, ∫QφdΩ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, ∫wVdΓ, ∫φMdΓ, L₂, L₂φ, L₂Q

using TimerOutputs, WriteVTK, XLSX
import Gmsh: gmsh

E = 10.92e6
ν = 0.3
h = 1e-0
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
# @timeit to "open msh file" gmsh.open("msh/patchtest_3.msh")
# @timeit to "get nodes" nodes_s = get𝑿ᵢ()

integrationOrder = 3
# ──────────────────────────────────────────────────────────
type_w = :tri3
type_φ = :tri3
type_Q = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type = eval(type_Q)
ndiv =8
ndiv_q = ndiv
XLSX.openxlsx("xls/square_33_8_8.xlsx", mode="w") do xf
 for ndiv_q = 8:32
 row = ndiv_q

@timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv_q.msh")
@timeit to "get nodes" nodes_q = get𝑿ᵢ()
@timeit to "get entities" entities = getPhysicalGroups()
xʷ = nodes_q.x
yʷ = nodes_q.y
zʷ = nodes_q.z
sp = RegularGrid(xʷ, yʷ, zʷ, n=3, γ=5)
nᵛ = length(nodes_q)
s = 1 / ndiv_q
s₁ = 1.5 * s * ones(nᵛ)
s₂ = 1.5 * s * ones(nᵛ)
s₃ = 1.5 * s * ones(nᵛ)
push!(nodes_q, :s₁ => s₁, :s₂ => s₂, :s₃ => s₃)

@timeit to "calculate error" begin
    @timeit to "get elements" elements_q = getElements(nodes_q, entities["Ω"], type, 10, sp)
    prescribe!(elements_q, :E => E, :ν => ν, :h => h, :Q₁ => Q₁, :Q₂ => Q₂)
    @timeit to "calculate shape functions" set𝝭!(elements_q)
end

@timeit to "open msh file" gmsh.open("msh/patchtest_$type_φ"*"_$ndiv.msh")
@timeit to "get nodes" nodes_φ = get𝑿ᵢ()
@timeit to "get entities" entities = getPhysicalGroups()

nᵠ = length(nodes_φ)
kᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)
kᵛᵛ = zeros(2 * nᵛ, 2 * nᵛ)
kᵛᵠ = zeros(2 * nᵛ, 2 * nᵠ)
fᵠ = zeros(2 * nᵠ)
fᵛ = zeros(2 * nᵛ)

@timeit to "calculate ∫κκdΩ" begin
    @timeit to "get elements" elements_Q = getElements(nodes_q, entities["Ω"], type, integrationOrder, sp)
    @timeit to "get elements" elements = getElements(nodes_φ, entities["Ω"], integrationOrder)
    prescribe!(elements, :E => E, :ν => ν, :h => h)
    prescribe!(elements_Q, :E => E, :ν => ν, :h => h)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    @timeit to "calculate shape functions" set𝝭!(elements_Q)
    𝑎ᵠᵠ = ∫κκdΩ => elements
    𝑎ᵛᵠ = ∫QφdΩ => (elements_Q, elements)
    𝑎ᵛᵛ = ∫QQdΩ => elements_Q
    # 𝑓ᵠ = ∫φmdΩ => elements
    @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)
    @timeit to "assemble" 𝑎ᵛᵛ(kᵛᵛ)
    @timeit to "assemble" 𝑎ᵛᵠ(kᵛᵠ)
    # @timeit to "assemble" 𝑓ᵠ(fᵠ)
end

@timeit to "calculate ∫αφφdΓ" begin
    @timeit to "get elements" elements_1 = getElements(nodes_φ, entities["Γ¹"], integrationOrder, normal=true)
    @timeit to "get elements" elements_2 = getElements(nodes_φ, entities["Γ²"], integrationOrder, normal=true)
    @timeit to "get elements" elements_3 = getElements(nodes_φ, entities["Γ³"], integrationOrder, normal=true)
    @timeit to "get elements" elements_4 = getElements(nodes_φ, entities["Γ⁴"], integrationOrder, normal=true)
    prescribe!(elements_1, :α => 1e8 * E, :g => w, :g₁ => φ₁, :g₂ => φ₂, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
    prescribe!(elements_2, :α => 1e8 * E, :g => w, :g₁ => φ₁, :g₂ => φ₂, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
    prescribe!(elements_3, :α => 1e8 * E, :g => w, :g₁ => φ₁, :g₂ => φ₂, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
    prescribe!(elements_4, :α => 1e8 * E, :g => w, :g₁ => φ₁, :g₂ => φ₂, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    @timeit to "calculate shape functions" set𝝭!(elements_4)
    𝑎ᵠ = ∫αφφdΓ => elements_1 ∪ elements_2 ∪ elements_3 ∪ elements_4
    @timeit to "assemble" 𝑎ᵠ(kᵠᵠ, fᵠ)
end

@timeit to "calculate error" begin
    @timeit to "get elements" elements_φ = getElements(nodes_φ, entities["Ω"], 10)
    prescribe!(elements_φ, :E => E, :ν => ν, :h => h, :φ₁ => φ₁, :φ₂ => φ₂)
    @timeit to "calculate shape functions" set𝝭!(elements_φ)
end

# ─── Defelection ──────────────────────────────────────────
@timeit to "open msh file" gmsh.open("msh/patchtest_$type_w"*"_$ndiv.msh")
@timeit to "get nodes" nodes_w = get𝑿ᵢ()
@timeit to "get entities" entities = getPhysicalGroups()

nʷ = length(nodes_w)
kʷʷ = zeros(nʷ, nʷ)
kᵛʷ = zeros(2 * nᵛ, nʷ)
kᵠʷ = zeros(2 * nᵠ, nʷ)
fʷ = zeros(nʷ)
@timeit to "calculate ∫Q∇wdΩ" begin
    @timeit to "get elements" elements = getElements(nodes_w, entities["Ω"], integrationOrder)
    prescribe!(elements, :E => E, :ν => ν, :h => h, :q => q)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    𝑎ᵛʷ = ∫Q∇wdΩ => (elements_Q, elements)
    𝑓ʷ = ∫wqdΩ => elements
    @timeit to "assemble" 𝑎ᵛʷ(kᵛʷ)
    @timeit to "assemble" 𝑓ʷ(fʷ)
end

@timeit to "calculate ∫αwwdΓ ∫QwdΓ" begin
    @timeit to "get elements" elements_1 = getElements(nodes_w, entities["Γ¹"], integrationOrder, normal=true)
    @timeit to "get elements" elements_2 = getElements(nodes_w, entities["Γ²"], integrationOrder, normal=true)
    @timeit to "get elements" elements_3 = getElements(nodes_w, entities["Γ³"], integrationOrder, normal=true)
    @timeit to "get elements" elements_4 = getElements(nodes_w, entities["Γ⁴"], integrationOrder, normal=true)
    @timeit to "get elements" elements_Q_1 = getElements(nodes_q, entities["Γ¹"], type, integrationOrder, sp, normal=true)
    @timeit to "get elements" elements_Q_2 = getElements(nodes_q, entities["Γ²"], type, integrationOrder, sp, normal=true)
    @timeit to "get elements" elements_Q_3 = getElements(nodes_q, entities["Γ³"], type, integrationOrder, sp, normal=true)
    @timeit to "get elements" elements_Q_4 = getElements(nodes_q, entities["Γ⁴"], type, integrationOrder, sp, normal=true)
    prescribe!(elements_1, :α => 1e8 * E, :g => w, :g₁ => φ₁, :g₂ => φ₂, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
    prescribe!(elements_2, :α => 1e8 * E, :g => w, :g₁ => φ₁, :g₂ => φ₂, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
    prescribe!(elements_3, :α => 1e8 * E, :g => w, :g₁ => φ₁, :g₂ => φ₂, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
    prescribe!(elements_4, :α => 1e8 * E, :g => w, :g₁ => φ₁, :g₂ => φ₂, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
    prescribe!(elements_Q_1, :g => w)
    prescribe!(elements_Q_2, :g => w)
    prescribe!(elements_Q_3, :g => w)
    prescribe!(elements_Q_4, :g => w)
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    @timeit to "calculate shape functions" set𝝭!(elements_4)
    @timeit to "calculate shape functions" set𝝭!(elements_Q_1)
    @timeit to "calculate shape functions" set𝝭!(elements_Q_2)
    @timeit to "calculate shape functions" set𝝭!(elements_Q_3)
    @timeit to "calculate shape functions" set𝝭!(elements_Q_4)
    𝑎ᵛ = ∫QwdΓ => (elements_Q_1 ∪ elements_Q_2 ∪ elements_Q_3 ∪ elements_Q_4, elements_1 ∪ elements_2 ∪ elements_3 ∪ elements_4)
    @timeit to "assemble" 𝑎ᵛ(kᵛʷ, fᵛ)
end

@timeit to "calculate error" begin
    @timeit to "get elements" elements_w = getElements(nodes_w, entities["Ω"], 10)
    prescribe!(elements_w, :E => E, :ν => ν, :h => h, :u => w)
    @timeit to "calculate shape functions" set𝝭!(elements_w)
end
# dᵠ = zeros(2*nᵠ)
# dᵛ = zeros(2*nᵛ)
# dʷ = zeros(nʷ)
# for node in nodes
#     x = node.x
#     y = node.y
#     z = node.z
#     dᵠ[2*node.𝐼-1] = φ₁(x,y,z)
#     dᵠ[2*node.𝐼]   = φ₂(x,y,z)
#     dᵛ[2*node.𝐼-1] = Q₁(x,y,z)
#     dᵛ[2*node.𝐼]   = Q₂(x,y,z)
#     dʷ[node.𝐼] = w(x,y,z)
# end
# println(kᵛᵛ*dᵛ)
# println(kᵛʷ*dʷ)
# println(kᵛᵛ*dᵛ + kᵛʷ*dʷ)
# println(kᵛʷ*ones(nʷ).-fᵛ)
# println(kᵠᵠ*dᵠ + kᵛᵠ'*dᵛ - fᵠ)
# println(kᵛᵛ*dᵛ + kᵛᵠ*dᵠ + kᵛʷ*dʷ - fᵛ)
# println(kᵛʷ'*dᵛ + kʷʷ*dʷ - fʷ)
# println(kᵛᵠ*dᵠ)
# println(kᵛʷ*dʷ)
# println(kᵛᵛ*dᵛ)
# println(kᵛᵛ*dᵛ + kᵛʷ*dʷ)

# println([kᵠᵠ kᵠʷ kᵛᵠ';kᵠʷ' kʷʷ kᵛʷ';kᵛᵠ kᵛʷ kᵛᵛ]*[dᵠ;dʷ;dᵛ] .- [fᵠ;fʷ;fᵛ])
@timeit to "solve" d = [kᵠᵠ kᵠʷ kᵛᵠ'; kᵠʷ' kʷʷ kᵛʷ'; kᵛᵠ kᵛʷ kᵛᵛ] \ [fᵠ; fʷ; fᵛ]
# println([kᵠᵠ kᵠʷ kᵛᵠ';kᵠʷ' kʷʷ kᵛʷ';kᵛᵠ kᵛʷ kᵛᵛ]*d .- [fᵠ;fʷ;fᵛ])
push!(nodes_φ, :d₁ => d[1:2:2*nᵠ], :d₂ => d[2:2:2*nᵠ])
push!(nodes_w, :d => d[2*nᵠ+1:2*nᵠ+nʷ])
push!(nodes_q, :q₁ => d[2*nᵠ+nʷ+1:2:end], :q₂ => d[2*nᵠ+nʷ+2:2:end])

@timeit to "calculate error" begin
    L₂_w = L₂(elements_w)
    L₂_φ = L₂φ(elements_φ)
    L₂_Q = L₂Q(elements_q)
end

points = zeros(3, nᵛ)
for node in nodes_q
    I = node.𝐼
    points[1, I] = node.x
    points[2, I] = node.y
    points[3, I] = node.z
end
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, [node.𝐼 for node in elm.𝓒]) for elm in elements]
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [node.𝐼 for node in elm.𝓒]) for elm in elements_q]
vtk_grid("vtk/square.vtu", points, cells) do vtk
    vtk["Q₁"] = [node.q₁ for node in nodes_q]
    vtk["Q₂"] = [node.q₂ for node in nodes_q]
    vtk["Q̄₁"] = [Q₁(node.x, node.y, node.z) for node in nodes_q]
    vtk["Q̄₂"] = [Q₂(node.x, node.y, node.z) for node in nodes_q]
end

# println(to)

println("L₂ error of w: ", L₂_w)
println("L₂ error of φ: ", L₂_φ)
println("L₂ error of Q: ", L₂_Q)
     # ──────────────────────────────────────────────────────────
     sheet = xf[1]
     XLSX.rename!(sheet, "new_sheet")
     sheet["A1"] = "type w"
     sheet["B1"] = "nʷ"
     sheet["C1"] = "type φ"
     sheet["D1"] = "nᵠ"
     sheet["E1"] = "type Q"
     sheet["F1"] = "nᵛ"
     sheet["G1"] = "L₂w"
     sheet["H1"] = "L₂φ"
     sheet["I1"] = "L₂Q"
     sheet["A$row"] = "$type_w"
     sheet["B$row"] = nʷ
     sheet["C$row"] = "$type_φ"
     sheet["D$row"] = nᵠ
     sheet["E$row"] = "$type_Q"
     sheet["F$row"] = nᵛ
     sheet["G$row"] = log10(L₂_w)
     sheet["H$row"] = log10(L₂_φ)
     sheet["I$row"] = log10(L₂_Q)
 end
    # end
    gmsh.finalize()
 end
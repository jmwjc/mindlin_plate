using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫wVdΓ, ∫φMdΓ, ∫αwwdΓ, ∫αφφdΓ, L₂, L₂φ

using WriteVTK
using TimerOutputs 
import Gmsh: gmsh

E = 1.0
ν = 0.3
h = 1e-1
Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))

w(x,y,z) = 1.0+x+y
w₁(x,y,z) = 1.0
w₂(x,y,z) = 1.0
w₁₁(x,y,z) = 0.0
w₂₂(x,y,z) = 0.0
φ₁(x,y,z) = 1.0+x+y
φ₂(x,y,z) = 1.0+x+y
φ₁₁(x,y,z)  = 1.0
φ₁₂(x,y,z)  = 1.0
φ₂₁(x,y,z)  = 1.0
φ₂₂(x,y,z)  = 1.0
φ₁₁₁(x,y,z)  = 0.0
φ₁₁₂(x,y,z)  = 0.0
φ₂₂₁(x,y,z)  = 0.0
φ₂₂₂(x,y,z)  = 0.0
φ₁₂₁(x,y,z)  = 0.0
φ₁₂₂(x,y,z)  = 0.0

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

n=16

gmsh.initialize()
@timeit to "open msh file" gmsh.open("msh/MorleysAcuteSkewPlate/MorleysAcuteSkewPlate_$n.msh")
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()

nʷ = length(nodes)
nᵠ = length(nodes)
kʷʷ = zeros(nʷ,nʷ)
kᵠᵠ = zeros(2*nᵠ,2*nᵠ)
kᵠʷ = zeros(2*nᵠ,nʷ)
fʷ = zeros(nʷ)
fᵠ = zeros(2*nᵠ)

@timeit to "calculate ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫wφdΩ" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"])
    prescribe!(elements, :E=>E, :ν=>ν, :h=>h, :q=>q, :m₁=>m₁, :m₂=>m₂)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    𝑎ʷʷ = ∫wwdΩ=>elements
    𝑎ᵠʷ = ∫φwdΩ=>elements
    𝑎ᵠᵠ = [
        ∫φφdΩ=>elements,
        ∫κκdΩ=>elements,
    ]
    𝑓ʷ = ∫wqdΩ=>elements
    𝑓ᵠ = ∫φmdΩ=>elements
    @timeit to "assemble" 𝑎ʷʷ(kʷʷ)
    @timeit to "assemble" 𝑎ᵠʷ(kᵠʷ)
    @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)
    @timeit to "assemble" 𝑓ʷ(fʷ)
    @timeit to "assemble" 𝑓ᵠ(fᵠ)

    global elements_domain = elements

end

@timeit to "calculate ∫αwwdΓ ∫αφφdΓ" begin
    @timeit to "get elements" elements_1 = getElements(nodes, entities["Γᵇ"])
    @timeit to "get elements" elements_2 = getElements(nodes, entities["Γʳ"])
    @timeit to "get elements" elements_3 = getElements(nodes, entities["Γᵗ"])
    @timeit to "get elements" elements_4 = getElements(nodes, entities["Γˡ"])

    prescribe!(elements_1, :α=>1e8*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_2, :α=>1e8*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_3, :α=>1e8*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_4, :α=>1e8*E, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)

    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    @timeit to "calculate shape functions" set𝝭!(elements_4)

    𝑎ʷ = ∫αwwdΓ=>elements_1∪elements_2∪elements_3∪elements_4
    𝑎ᵠ = ∫αφφdΓ=>elements_1∪elements_2∪elements_3∪elements_4
    @timeit to "assemble" 𝑎ʷ(kʷʷ, fʷ)
    @timeit to "assemble" 𝑎ᵠ(kᵠᵠ, fᵠ)
end

@timeit to "solve" d = [kᵠᵠ kᵠʷ;kᵠʷ' kʷʷ]\[fᵠ;fʷ]
push!(nodes,:d=>d[2*nᵠ+1:end], :d₁=>d[1:2:2*nᵠ], :d₂=>d[2:2:2*nᵠ])

@timeit to "calculate error" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"], 10)
    prescribe!(elements, :E=>E, :ν=>ν, :h=>h, :u=>w, :φ₁=>φ₁, :φ₂=>φ₂)
    @timeit to "calculate shape functions" set𝝭!(elements)
    L₂_w = L₂(elements)
    L₂_φ = L₂φ(elements)
end
 
gmsh.finalize()

println(to)

println("L₂ error of w: ", L₂_w)
println("L₂ error of φ: ", L₂_φ)


# 图--------------------------------------------------------------------------------

# 坐标------------------------------------------------------------------------------
# nₚ = length(nodes)
# points = zeros(3,nₚ)
# for (i,node) in enumerate(nodes)
#     points[1,i] = node.x
#     points[2,i] = node.y
#     points[3,i] = node.d*4
#     # points[3,i] = us[i]*4
# end

nₚ = length(nodes)
points = zeros(3,nₚ)
for (i,node) in enumerate(nodes)
    points[1,i] = node.x
    points[2,i] = node.y
    points[3,i] = node.d
    # points[3,i] = us[i]*4
end

# 二维------------------------------------------------------------------------------
# xs = [node.x for node in nodes]'
# ys = [node.y for node in nodes]'
# zs = [node.z for node in nodes]'
# points = [xs; ys; zs]
# cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP,[x.𝐼 for x in elm.𝓒]) for elm in elements["Ω"]]
# # vtk_grid("./vtk/hmd_2d/error/non_uniform_Tri3_"*string(ndiv)*".vtu",points,cells) do vtk
# vtk_grid("./vtk/hmd_2d/Tri3_d_"*string(ndiv)*".vtu",points,cells) do vtk
#     vtk["d"] = [node.d for node in nodes]
#     # vtk["精确解"] = us
# end

# fₓ,fₜ,fₓₓ,fₜₜ = truncation_error(elements["Ω"],nₚ)
# println(fₓ)
# println(fₜ)
# println(fₛ)

# 三维------------------------------------------------------------------------------

# # cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE_STRIP, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ω"]]
# # vtk_grid("./vtk/circular_tri3_"*string(ndiv), points, cells) do vtk
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements_domain]


# 三维误差-------------------------------------------------------------------------------------------------
# vtk_grid("./vtk/circular_tri3_" * string(ndiv) * ".vtu", points, cells;
#          ascii=true, append=false, compress=false) do vtk


#     vtk["L2_w", WriteVTK.VTKFieldData()] = [L₂_w]
#     vtk["L2_phi", WriteVTK.VTKFieldData()] = [L₂_φ]
# end
# -------------------------------------------------------------------------------------------------

# 三维变形------------------------------------------------------------------------------
# vtk_grid("./vtk/circular_Clamped_tri3_" * string(ndiv) * ".vtu", points, cells;
#          ascii=true, append=false, compress=false) do vtk
vtk_grid("./vtk/MorleysAcuteSkewPlate_tri3_$n.vtu", points, cells;
         ascii=true, append=false, compress=false) do vtk

    # 挠度 w
    vtk["w"] = [node.d for node in nodes]
    # 转角 φ₁
    vtk["φ₁"] = [node.d₁ for node in nodes]
    # 转角 φ₂  
    vtk["φ₂"] = [node.d₂ for node in nodes]
end
# -------------------------------------------------------------------------------------------------

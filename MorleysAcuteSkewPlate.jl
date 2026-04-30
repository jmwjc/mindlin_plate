using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫wVdΓ, ∫φMdΓ, ∫αwwdΓ, ∫αφφdΓ, L₂, L₂φ

using WriteVTK
using TimerOutputs 
import Gmsh: gmsh

# Razzaque 60° skew plate (Katili 1993, Fig.19 / Table VIII)
L = 100.0
E = 1085.0
ν = 0.31
h = 0.001 * L            # = 0.1
Dᵇ = E*h^3/12/(1-ν^2)    # 后面做 w_c 无量纲化要用到
Dˢ = 5/6*E*h/(2*(1+ν))   # 仅作记录：算子内部也固定用 κ=5/6

# uniform transverse load
fz = 1.0
q(x,y,z)  = fz
m₁(x,y,z) = 0.0
m₂(x,y,z) = 0.0

# hard targets on Γᵇ and Γᵗ (先用这些名字保证后续 prescribe! 不报错)
w(x,y,z)  = 0.0
φ₁(x,y,z) = 0.0
φ₂(x,y,z) = 0.0

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
    @timeit to "get elements" elements_b = getElements(nodes, entities["Γᵇ"])
    @timeit to "get elements" elements_t = getElements(nodes, entities["Γᵗ"])

    α = 1e8 * E
    prescribe!(elements_b, :α=>α, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_t, :α=>α, :g=>w, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)

    @timeit to "calculate shape functions" set𝝭!(elements_b)
    @timeit to "calculate shape functions" set𝝭!(elements_t)

    elements_Γ = elements_b ∪ elements_t

    𝑎ʷ = ∫αwwdΓ=>elements_Γ
    𝑎ᵠ = ∫αφφdΓ=>elements_Γ

    @timeit to "assemble" 𝑎ʷ(kʷʷ, fʷ)
    @timeit to "assemble" 𝑎ᵠ(kᵠᵠ, fᵠ)
end

@timeit to "solve" d = [kᵠᵠ kᵠʷ;kᵠʷ' kʷʷ]\[fᵠ;fʷ]
push!(nodes,:d=>d[2*nᵠ+1:end], :d₁=>d[1:2:2*nᵠ], :d₂=>d[2:2:2*nᵠ])

@timeit to "postprocess center wc" begin
    # Physical Point "𝐴" 在该 msh 里是板中心点
    elements_A = getElements(nodes, entities["𝐴"])
    I_A = elements_A[1].𝓒[1].𝐼

    w_center = nodes[I_A].d
    w_c = w_center * 1e3 * Dᵇ / (fz * L^4)   # Table VIII 定义

    println("w(center) = ", w_center)
    println("w_c = ", w_c, "   (ref: 0.7945)")
end

gmsh.finalize()

println(to)


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
    points[3,i] = node.d/20
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

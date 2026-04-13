using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, L₂, L₂φ, L₂Q, ∫wVdΓ, ∫φMdΓ, ∫QφdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ
using WriteVTK
using LinearAlgebra
using TimerOutputs
import Gmsh: gmsh

# -----------------------------
# case switch (reference solution from 圆板精確解.md)
# -----------------------------

    E = 10.92
    ν = 0.3
    h = 1e0
    R = 5.0
    fz = 1.0
    α = 1e8 * E 
    Dᵇ = E*h^3/12/(1-ν^2)
    Dˢ = 5/6*E*h/(2*(1+ν))
    
 



# 解析解------------------------------------------------------------------------------

r(x, y) = sqrt(x^2 + y^2)
ξ(x, y) = r(x, y) / R
cosθ(x, y) = x/r(x,y)
sinθ(x, y) = y/r(x,y)

w(x, y, z) = fz * R^4 / (64 * Dᵇ) * (1 - ξ(x, y)^2) * ((1 - ξ(x, y)^2) + 8 * h^2 / (3 * (5/6) * R^2 * (1 - ν)))


φᵣ(x, y, z) = fz * r(x, y) * (r(x, y)^2 - R^2) / (16 * Dᵇ)
φ₁(x,y,z) = φᵣ(x,y,z)*x/r(x,y)
φ₂(x,y,z) = φᵣ(x,y,z)*y/r(x,y)
φ₁₁(x,y,z) = φᵣ(x,y,z)*cos²θ - φᵣ(x,y,z)*sin²θ
φ₁₂(x,y,z) = φᵣ(x,y,z)*sinθ*cosθ
φ₂₂(x,y,z) = φᵣ(x,y,z)*sin²θ - φᵣ(x,y,z)*cos²θ
φ₁₁₁(x,y,z) = φᵣ(x,y,z)*cos³θ - 3*φᵣ(x,y,z)*sin²θ*cosθ
φ₁₁₂(x,y,z) = φᵣ(x,y,z)*sinθ*cos²θ - φᵣ(x,y,z)*sin³θ
φ₂₂₁(x,y,z) = φᵣ(x,y,z)*sin²θ*cosθ - φᵣ(x,y,z)*cos³θ
φ₂₂₂(x,y,z) = φᵣ(x,y,z)*sin³θ - 3*φᵣ(x,y,z)*sinθ*cos²θ
φ₁₂₁(x,y,z) = φᵣ(x,y,z)*sin²θ*cosθ + φᵣ(x,y,z)*cos³θ
φ₁₂₂(x,y,z) = φᵣ(x,y,z)*sin³θ + φᵣ(x,y,z)*sinθ*cos²θ

κᵣᵣ(x, y, z) = -fz * R^2 * (3 * ξ(x, y)^2 - 1) / (16 * Dᵇ)
κᵩᵩ(x, y, z) = -fz * R^2 * (ξ(x, y)^2 - 1) / (16 * Dᵇ)

γᵣ(x, y, z) = -fz * r(x, y) / Dˢ

Mᵣᵣ(x, y, z) = -fz * R^2 / 16 * ((3 + ν) * ξ(x, y)^2 - (1 + ν))
Mᵩᵩ(x, y, z) = -fz * R^2 / 16 * ((1 + 3ν) * ξ(x, y)^2 - (1 + ν))

M₁₁(x,y,z) = Mᵣᵣ(x,y,z) * cos²θ + Mᵩᵩ(x,y,z) * sin²θ
M₂₂(x,y,z) = Mᵣᵣ(x,y,z) * sin²θ + Mᵩᵩ(x,y,z) * cos²θ
M₁₂(x,y,z) = (Mᵣᵣ(x,y,z) - Mᵩᵩ(x,y,z)) * sinθ * cosθ
M₁₁₁(x,y,z) = fz * x * (1 - ν) / 2 * (1 - ξ(x,y)^2)
M₁₂₂(x,y,z) = fz * x * y * (1 - ν) / (2 * R^2) * ξ(x,y)
M₁₂₁(x,y,z) = fz * x * y * (1 - ν) / (2 * R^2) * ξ(x,y)
M₂₂₂(x,y,z) = fz * y * (1 - ν) / 2 * (1 - ξ(x,y)^2)
m₁(x,y,z) = M₁₁₁(x,y,z) + M₁₂₂(x,y,z) - Q₁(x,y,z)
m₂(x,y,z) = M₁₂₁(x,y,z) + M₂₂₂(x,y,z) - Q₂(x,y,z)

Qᵣ(x, y, z) = -fz * r(x, y) / 2
Q₁(x, y, z) = -fz * x / 2
Q₂(x, y, z) = -fz * y / 2
Q₁₁(x,y,z) = -fz / 2  
Q₂₂(x,y,z) = -fz / 2  

# --------------------------------------------------------------------------------


# 自然
# （示意图边界验收时用零边界；解析解对比时也依然是 clamped 零边界）
# w0(x, y, z) = 0.0
# φ10(x, y, z) = 0.0
# φ20(x, y, z) = 0.0

# V(x, y, z) = 0.0
# M₁(x, y, z) = 0.0
# M₂(x, y, z) = 0.0



const to = TimerOutput()
n = 15
gmsh.initialize()
# @timeit to "open msh file" gmsh.open("./msh/circular.msh")
# @timeit to "open msh file" gmsh.open("./msh/circular_tri3_9.msh")
@timeit to "open msh file" gmsh.open("./msh/circular_tri3_$n.msh")
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()

nʷ = length(nodes)
nᵠ = length(nodes)
nᵛ = length(nodes)

kʷʷ = zeros(nʷ, nʷ)
kᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)
kᵛᵛ = zeros(2*nᵛ,2*nᵛ)
kᵠʷ = zeros(2 * nᵠ, nʷ)
kᵛʷ = zeros(2*nᵛ,nʷ)
kᵛᵠ = zeros(2*nᵛ,2*nᵠ)

fʷ = zeros(nʷ)
fᵠ = zeros(2 * nᵠ)
fᵛ = zeros(2*nᵛ)

integrationOrder = 3 
@timeit to "assemble domain" begin
    elements = getElements(nodes, entities["Ω"], integrationOrder)
    elements_Γᵇ = getElements(nodes, entities["Γᵇ"], integrationOrder, normal=true)
    elements_Γᵉ = getElements(nodes, entities["Γᵉ"], integrationOrder, normal=true)
    elements_Γˡ = getElements(nodes, entities["Γˡ"], integrationOrder, normal=true)
    elements_Γ = elements_Γᵇ ∪ elements_Γᵉ ∪ elements_Γˡ
    prescribe!(elements, :E => E, :ν => ν, :h => h, :q => fz)
    set∇𝝭!(elements)
     set𝝭!(elements_Γᵇ)
    set𝝭!(elements_Γᵉ)
    set𝝭!(elements_Γˡ)

    # 𝑎ʷʷ = ∫wwdΩ => elements
    # 𝑎ᵠʷ = ∫φwdΩ => elements
    # 𝑎ᵠᵠ = [
    #     ∫φφdΩ => elements,
    #     ∫κκdΩ => elements,
    # ]
    # 𝑓ʷ = ∫wqdΩ => elements
    # # 𝑓ᵠ = ∫φmdΩ => elements
    # @timeit to "assemble" 𝑎ʷʷ(kʷʷ)
    # @timeit to "assemble" 𝑎ᵠʷ(kᵠʷ)
    # @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)
    # @timeit to "assemble" 𝑓ʷ(fʷ)
    # # @timeit to "assemble" 𝑓ᵠ(fᵠ)

    𝑎ᵠᵠ = ∫κκdΩ=>elements
    𝑎ᵛᵠ = ∫QφdΩ=>(elements, elements)
    𝑎ᵛᵛ = ∫QQdΩ=>elements
    𝑎ᵛʷ = [
        ∫∇QwdΩ=>(elements, elements),
        ∫QwdΓ=>(elements_Γ, elements_Γ),
    ]
    𝑓ʷ = ∫wqdΩ => elements

    @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)
    @timeit to "assemble" 𝑎ᵛᵠ(kᵛᵠ)
    @timeit to "assemble" 𝑎ᵛᵛ(kᵛᵛ)
    @timeit to "assemble" 𝑎ᵛʷ(kᵛʷ)
    @timeit to "assemble" 𝑓ʷ(fʷ)

    global elements_domain = elements 

end



@timeit to "Γᵉ" begin
    elements = getElements(nodes, entities["Γᵉ"])
    set𝝭!(elements)

    prescribe!(elements, 
        :α => α,
        :g => w,
        :g₁ => φ₁,
        :g₂ => φ₂,
        :n₁₁ => 1.0,
        :n₁₂ => 0.0,
        :n₂₂ => 1.0,
    )
    𝑎ʷ = ∫αwwdΓ => elements
    𝑎ᵠ = ∫αφφdΓ => elements
    @timeit to "assemble" 𝑎ʷ(kʷʷ,fʷ)
    @timeit to "assemble" 𝑎ᵠ(kᵠᵠ,fᵠ)
 
end

@timeit to "Γˡ" begin
    elements = getElements(nodes, entities["Γˡ"],normal=true)
    set𝝭!(elements)

    prescribe!(elements, 
        :α => α,
        :g₁ => φ₁,
        :g₂ => φ₂,
        :n₁₁ => 1.0,
        :n₁₂ => 0.0,
        :n₂₂ => 0.0,
    )
    𝑎ᵠ = ∫αφφdΓ => elements
    @timeit to "assemble" 𝑎ᵠ(kᵠᵠ,fᵠ)
end

@timeit to "Γᵇ" begin
    elements = getElements(nodes, entities["Γᵇ"],normal=true)
    set𝝭!(elements)

    prescribe!(elements, 
        :α => α,
        :g₁ => φ₁,
        :g₂ => φ₂,
        :n₁₁ => 0.0,
        :n₁₂ => 0.0,
        :n₂₂ => 1.0,
    )
    𝑎ᵠ = ∫αφφdΓ => elements
    @timeit to "assemble" 𝑎ᵠ(kᵠᵠ,fᵠ)
end



@timeit to "solve" d = [kᵠᵠ kᵠʷ kᵛᵠ'; kᵠʷ' kʷʷ kᵛʷ'; kᵛᵠ kᵛʷ kᵛᵛ] \ [fᵠ; fʷ; fᵛ]

push!(nodes, 
    :d => d[2*nᵠ+1:2*nᵠ+nʷ],       # w
    :d₁ => d[1:2:2*nᵠ],            # φ₁
    :d₂ => d[2:2:2*nᵠ],            # φ₂
    :q₁ => d[2*nᵠ+nʷ+1:2:end],     # Q₁
    :q₂ => d[2*nᵠ+nʷ+2:2:end]      # Q₂
)

@timeit to "calculate error" begin
    elements_err = getElements(nodes, entities["Ω"], 10)
    prescribe!(elements_err, :E => E, :ν => ν, :h => h, :u => w, :φ₁ => φ₁, :φ₂ => φ₂, :Q₁ => Q₁, :Q₂ => Q₂)
    set𝝭!(elements_err)
    global L₂_w = L₂(elements_err)
    global L₂_φ = L₂φ(elements_err)
    global L₂_Q = L₂Q(elements_err)
end


gmsh.finalize()

println(to)

println("α penalty: ", α)
println("L₂ error of w: ", L₂_w)
println("L₂ error of φ: ", L₂_φ)
println("L₂ error of Q: ", L₂_Q)
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
    points[3,i] = node.d/18
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
vtk_grid("./vtk/circular_Clamped_mix_tri3_$n.vtu", points, cells;
         ascii=true, append=false, compress=false) do vtk

    # 挠度 w
    vtk["w"] = [node.d for node in nodes]
    # 转角 φ₁
    vtk["phi_1"] = [node.d₁ for node in nodes]
    # 转角 φ₂  
    vtk["phi_2"] = [node.d₂ for node in nodes]

    vtk["Q_1"] = [node.q₁ for node in nodes]
      
    vtk["Q_2"] = [node.q₂ for node in nodes]
end
# -------------------------------------------------------------------------------------------------
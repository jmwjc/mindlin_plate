using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, L₂, L₂φ, L₂Q, ∫wVdΓ, ∫φMdΓ, ∫QφdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ
using WriteVTK
using LinearAlgebra
using TimerOutputs
using XLSX
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
φ₁(x,y,z) = r(x,y) ≈ 0.0 ? 0.0 : φᵣ(x,y,z)*x/r(x,y)
φ₂(x,y,z) = r(x,y) ≈ 0.0 ? 0.0 : φᵣ(x,y,z)*y/r(x,y)
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


const to = TimerOutput()
gmsh.initialize()
integrationOrder = 5
# ──────────────────────────────────────────────────────────
type_w = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_φ = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_Q = :tri3
type_M = :(PiecewisePolynomial{:Linear2D})
ndiv_φ = 16
ndiv_w = 16
ndiv = ndiv_φ
XLSX.openxlsx("xls/circular_clamped_$(ndiv_φ)_tri3_$(ndiv_w).xlsx", mode="w") do xf
    # for ndiv = ndiv_w:25
        row = ndiv      
# ─── Deflection W ─────────────────────────────────────────
@timeit to "open msh file" gmsh.open("./msh/circular_tri3_$ndiv_w.msh")
@timeit to "get nodes" nodes_w = get𝑿ᵢ()
xʷ = nodes_w.x
yʷ = nodes_w.y
zʷ = nodes_w.z
sp_w = RegularGrid(xʷ,yʷ,zʷ,n = 3,γ = 5)
nʷ = length(nodes_w)
s = R/ndiv_w
s₁ = 1.5*s*ones(nʷ)
s₂ = 1.5*s*ones(nʷ)
s₃ = 1.5*s*ones(nʷ)
push!(nodes_w,:s₁=>s₁,:s₂=>s₂,:s₃=>s₃)
# ─── Rotation Φ ───────────────────────────────────────────
@timeit to "open msh file" gmsh.open("./msh/circular_tri3_$ndiv_φ.msh")
@timeit to "get nodes" nodes_φ = get𝑿ᵢ()
xᵠ = nodes_φ.x
yᵠ = nodes_φ.y
zᵠ = nodes_φ.z
sp_φ = RegularGrid(xᵠ,yᵠ,zᵠ,n = 3,γ = 5)
nᵠ = length(nodes_φ)
s = R/ndiv_φ
s₁ = 1.5*s*ones(nᵠ)
s₂ = 1.5*s*ones(nᵠ)
s₃ = 1.5*s*ones(nᵠ)
push!(nodes_φ,:s₁=>s₁,:s₂=>s₂,:s₃=>s₃)
# ─── Shear ────────────────────────────────────────────────
@timeit to "open msh file" gmsh.open("./msh/circular_tri3_$ndiv.msh")
@timeit to "get nodes" nodes = get𝑿ᵢ()
@timeit to "get entities" entities = getPhysicalGroups()

nˢ = length(nodes)
kᵠᵠ = zeros(2*nᵠ,2*nᵠ)
kʷʷ = zeros(nʷ,nʷ)
kˢˢ = zeros(2*nˢ,2*nˢ)
kˢᵠ = zeros(2*nˢ,2*nᵠ)
kˢʷ = zeros(2*nˢ,nʷ)
kᵠʷ = zeros(2*nᵠ,nʷ)
fᵠ = zeros(2*nᵠ)
fʷ = zeros(nʷ)
fˢ = zeros(2*nˢ)

    # ── Domain elements (Q/FEM + w/mf + φ/mf) ─────────────────────────────
    elements_q = getElements(nodes,   entities["Ω"], integrationOrder)
    elements_w = getElements(nodes_w, entities["Ω"], eval(type_w), integrationOrder, sp_w)
    elements_φ = getElements(nodes_φ, entities["Ω"], eval(type_φ), integrationOrder, sp_φ)

    prescribe!(elements_q, :E => E, :ν => ν, :h => h)
    prescribe!(elements_w, :E => E, :ν => ν, :h => h, :q => fz)
    prescribe!(elements_φ, :E => E, :ν => ν, :h => h)

    set∇𝝭!(elements_q)
    set𝝭!(elements_w)
    set∇𝝭!(elements_φ)

    # ── Boundary elements for ∫QwdΓ (Route B: Γ = Γᵇ ∪ Γᵉ ∪ Γˡ) ─────────
    elements_q_Γᵇ = getElements(nodes, entities["Γᵇ"], integrationOrder, normal=true)
    elements_q_Γᵉ = getElements(nodes, entities["Γᵉ"], integrationOrder, normal=true)
    elements_q_Γˡ = getElements(nodes, entities["Γˡ"], integrationOrder, normal=true)

    elements_w_Γᵇ = getElements(nodes_w, entities["Γᵇ"], eval(type_w), integrationOrder, sp_w, normal=true)
    elements_w_Γᵉ = getElements(nodes_w, entities["Γᵉ"], eval(type_w), integrationOrder, sp_w, normal=true)
    elements_w_Γˡ = getElements(nodes_w, entities["Γˡ"], eval(type_w), integrationOrder, sp_w, normal=true)

    set𝝭!(elements_q_Γᵇ); set𝝭!(elements_q_Γᵉ); set𝝭!(elements_q_Γˡ)
    set𝝭!(elements_w_Γᵇ); set𝝭!(elements_w_Γᵉ); set𝝭!(elements_w_Γˡ)

    elements_q_Γ = elements_q_Γᵇ ∪ elements_q_Γᵉ ∪ elements_q_Γˡ
    elements_w_Γ = elements_w_Γᵇ ∪ elements_w_Γᵉ ∪ elements_w_Γˡ

    # ── Operators (3-field mix, shear named ˢ) ───────────────────────────
    𝑎ᵠᵠ = ∫κκdΩ => elements_φ
    𝑎ˢᵠ = ∫QφdΩ => (elements_q, elements_φ)
    𝑎ˢˢ = ∫QQdΩ => elements_q
    𝑎ˢʷ = [
        ∫∇QwdΩ => (elements_q, elements_w),
        ∫QwdΓ  => (elements_q_Γ, elements_w_Γ),
    ]
    𝑓ʷ = ∫wqdΩ => elements_w

    @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)
    @timeit to "assemble" 𝑎ˢᵠ(kˢᵠ)
    @timeit to "assemble" 𝑎ˢˢ(kˢˢ)
    @timeit to "assemble" 𝑎ˢʷ(kˢʷ)
    @timeit to "assemble" 𝑓ʷ(fʷ) 

nₑ = length(elements_q)
nᵐ = nₑ * ApproxOperator.get𝑛𝑝(eval(type_M)(𝑿ᵢ[], 𝑿ₛ[]))

@timeit to "Γᵉ" begin
elements_w = getElements(nodes_w, entities["Γᵉ"], eval(type_w), integrationOrder, sp_w, normal=true)
set𝝭!(elements_w)

prescribe!(elements_w,
:α => α,
:g => w,
)
𝑎ʷ = ∫αwwdΓ => elements_w
@timeit to "assemble" 𝑎ʷ(kʷʷ, fʷ)

elements_φ = getElements(nodes_φ, entities["Γᵉ"], eval(type_φ), integrationOrder, sp_φ, normal=true)
set𝝭!(elements_φ)

prescribe!(elements_φ,
:α => α,
:g₁ => φ₁,
:g₂ => φ₂,
:n₁₁ => 1.0,
:n₁₂ => 0.0,
:n₂₂ => 1.0,
)
𝑎ᵠ = ∫αφφdΓ => elements_φ
@timeit to "assemble" 𝑎ᵠ(kᵠᵠ, fᵠ)
end

@timeit to "Γˡ" begin
elements = getElements(nodes_φ, entities["Γˡ"], eval(type_φ), integrationOrder, sp_φ, normal=true)
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
@timeit to "assemble" 𝑎ᵠ(kᵠᵠ, fᵠ)
end

@timeit to "Γᵇ" begin
elements = getElements(nodes_φ, entities["Γᵇ"], eval(type_φ), integrationOrder, sp_φ, normal=true)
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
@timeit to "assemble" 𝑎ᵠ(kᵠᵠ, fᵠ)
end



@timeit to "solve" d = [kᵠᵠ kᵠʷ kˢᵠ'; kᵠʷ' kʷʷ kˢʷ'; kˢᵠ kˢʷ kˢˢ] \ [fᵠ; fʷ; fˢ]

push!(nodes_φ,
    :d₁ => d[1:2:2*nᵠ],
    :d₂ => d[2:2:2*nᵠ],
)

push!(nodes_w,
    :d => d[2*nᵠ+1:2*nᵠ+nʷ],
)

push!(nodes,
    :q₁ => d[2*nᵠ+nʷ+1:2:2*nᵠ+nʷ+2*nˢ],
    :q₂ => d[2*nᵠ+nʷ+2:2:2*nᵠ+nʷ+2*nˢ],
)

@timeit to "calculate error" begin
    elements_w_err = getElements(nodes_w, entities["Ω"], eval(type_w), 10, sp_w)
    prescribe!(elements_w_err, :E => E, :ν => ν, :h => h, :u => w)
    set𝝭!(elements_w_err)

    elements_φ_err = getElements(nodes_φ, entities["Ω"], eval(type_φ), 10, sp_φ)
    prescribe!(elements_φ_err, :E => E, :ν => ν, :h => h, :φ₁ => φ₁, :φ₂ => φ₂)
    set𝝭!(elements_φ_err)

    elements_q_err = getElements(nodes, entities["Ω"], 10)
    prescribe!(elements_q_err, :E => E, :ν => ν, :h => h, :Q₁ => Q₁, :Q₂ => Q₂)
    set𝝭!(elements_q_err)

    global L₂_w = L₂(elements_w_err)
    global L₂_φ = L₂φ(elements_φ_err)
    global L₂_Q = L₂Q(elements_q_err)
end




println(to)

println("α penalty: ", α)
println("nʷ=$nʷ, nᵠ=$nᵠ, nˢ=$nˢ, nᵐ=$nᵐ")
println("L₂ error of w: ", L₂_w)
println("L₂ error of φ: ", L₂_φ)
println("L₂ error of Q: ", L₂_Q)

wˢ = zeros(nˢ)
φ₁ˢ = zeros(nˢ)
φ₂ˢ = zeros(nˢ)
for xᵢ in nodes
    indices = sp_w(xᵢ.x,xᵢ.y,xᵢ.z)
    ξᵢ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=length(indices)), Dict([:x=>(1,[xᵢ.x]),:y=>(1,[xᵢ.y]),:z=>(1,[xᵢ.z])]))
    𝓒 = nodes_w[[indices...]]
    𝓖 = [ξᵢ]
    element = [eval(type_w)(𝓒,𝓖)]
    set𝝭!(element)
    wᵢ = 0.0
    N = ξᵢ[:𝝭]
    for (i,xⱼ) in enumerate(𝓒)
        wᵢ += N[i]*xⱼ.d
    end
    wˢ[xᵢ.𝐼] = wᵢ

    indices = sp_φ(xᵢ.x,xᵢ.y,xᵢ.z)
    ξᵢ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=length(indices)), Dict([:x=>(1,[xᵢ.x]),:y=>(1,[xᵢ.y]),:z=>(1,[xᵢ.z])]))
    𝓒 = nodes_φ[[indices...]]
    𝓖 = [ξᵢ]
    element = [eval(type_φ)(𝓒,𝓖)]
    set𝝭!(element)
    φ₁ᵢ = 0.0
    φ₂ᵢ = 0.0
    N = ξᵢ[:𝝭]
    for (i,xⱼ) in enumerate(𝓒)
        φ₁ᵢ += N[i]*xⱼ.d₁
        φ₂ᵢ += N[i]*xⱼ.d₂
    end
    φ₁ˢ[xᵢ.𝐼] = φ₁ᵢ
    φ₂ˢ[xᵢ.𝐼] = φ₂ᵢ
end
push!(nodes,:w=>wˢ,:φ₁=>φ₁ˢ,:φ₂=>φ₂ˢ)

nₚ = length(nodes)
points = zeros(3,nₚ)
for (i,node) in enumerate(nodes)
    points[1,i] = node.x
    points[2,i] = node.y
    points[3,i] = node.w/15
    # points[3,i] = 0.0
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
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, [xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements_q]


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
vtk_grid("./vtk/circular_clamped_$(ndiv_φ)_tri3_$(ndiv_w).vtu", points, cells;
         ascii=true, append=false, compress=false) do vtk

   
    vtk["w"] = [node.w for node in nodes]
    vtk["φ₁"] = [node.φ₁ for node in nodes]
    vtk["φ₂"] = [node.φ₂ for node in nodes]
    vtk["q₁"] = [node.q₁ for node in nodes]
    vtk["q₂"] = [node.q₂ for node in nodes]
    vtk["qᵣ"] = [r(node.x,node.y) ≈ 0.0 ? 0.0 : node.q₁*node.x/r(node.x,node.y) + node.q₂*node.y/r(node.x,node.y) for node in nodes]

    vtk["w̄"] = [w(node.x,node.y,node.z) for node in nodes]
    vtk["φ̄₁"] = [φ₁(node.x,node.y,node.z) for node in nodes]
    vtk["φ̄₂"] = [φ₂(node.x,node.y,node.z) for node in nodes]
    vtk["q̄₁"] = [Q₁(node.x,node.y,node.z) for node in nodes]
    vtk["q̄₂"] = [Q₂(node.x,node.y,node.z) for node in nodes]
    vtk["q̄ᵣ"] = [Qᵣ(node.x,node.y,node.z) for node in nodes]
end

# ──────────────────────────────────────────────────────────
     sheet = xf[1]
     XLSX.rename!(sheet, "new_sheet")
     sheet["A1"] = "type w"
     sheet["B1"] = "nʷ"
     sheet["C1"] = "type φ"
     sheet["D1"] = "nᵠ"
     sheet["E1"] = "type Q"
     sheet["F1"] = "nˢ"
     sheet["G1"] = "type M"
     sheet["H1"] = "nᵐ"
     sheet["I1"] = "L₂w"
     sheet["J1"] = "L₂φ"
     sheet["K1"] = "L₂Q"
     sheet["A$row"] = "$type_w"
     sheet["B$row"] = nʷ
     sheet["C$row"] = "$type_φ"
     sheet["D$row"] = nᵠ
     sheet["E$row"] = "$type_Q"
     sheet["F$row"] = nˢ
     sheet["G$row"] = "$type_M"
     sheet["H$row"] = nᵐ
     sheet["I$row"] = log10(L₂_w)
     sheet["J$row"] = log10(L₂_φ)
     sheet["K$row"] = log10(L₂_Q)

    # end
    end
    gmsh.finalize()
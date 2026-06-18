using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements, getPiecewiseElements, getPiecewiseBoundaryElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫QQdΩ, ∫∇QwdΩ, ∫QwdΓ, ∫QφdΩ, ∫MMdΩ, ∫∇MφdΩ, ∫MφdΓ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, ∫wwdΩ , ∫wVdΓ, ∫φMdΓ

using TimerOutputs, LinearAlgebra, XLSX
import Gmsh: gmsh
include("cal_area_support.jl")


E = 1e6
ν = 0.3
h = 1e-3
Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))
integrationOrder = 2

# ──────────────────────────────────────────────────────────
w(x, y, z) = 1.0 + x + y
w₁(x, y, z) = 1.0
w₂(x, y, z) = 1.0
w₁₁(x, y, z) = 0.0
w₂₂(x, y, z) = 0.0
φ₁(x, y, z) = 1.0 + x + y
φ₂(x, y, z) = 1.0 + x + y
φ₁₁(x, y, z) = 1.0
φ₁₂(x, y, z) = 1.0
φ₂₁(x, y, z) = 1.0
φ₂₂(x, y, z) = 1.0
φ₁₁₁(x, y, z) = 0.0
φ₁₁₂(x, y, z) = 0.0
φ₂₂₁(x, y, z) = 0.0
φ₂₂₂(x, y, z) = 0.0
φ₁₂₁(x, y, z) = 0.0
φ₁₂₂(x, y, z) = 0.0

# r = 1
# w(x,y,z) = (x+y)^r
# w₁(x,y,z) = r*(x+y)^abs(r-1)
# w₂(x,y,z) = r*(x+y)^abs(r-1)
# w₁₁(x,y,z) = r*(r-1)*(x+y)^abs(r-2)
# w₂₂(x,y,z) = r*(r-1)*(x+y)^abs(r-2)
# φ₁(x,y,z) = r*(x+y)^abs(r-1)
# φ₂(x,y,z) = r*(x+y)^abs(r-1)
# φ₁₁(x,y,z)  = r*(r-1)*(x+y)^abs(r-2)
# φ₁₂(x,y,z)  = r*(r-1)*(x+y)^abs(r-2)
# φ₂₁(x,y,z)  = r*(r-1)*(x+y)^abs(r-2)
# φ₂₂(x,y,z)  = r*(r-1)*(x+y)^abs(r-2)
# φ₁₁₁(x,y,z)  = r*(r-1)*(r-2)*(x+y)^abs(r-3)
# φ₁₁₂(x,y,z)  = r*(r-1)*(r-2)*(x+y)^abs(r-3)
# φ₂₂₁(x,y,z)  = r*(r-1)*(r-2)*(x+y)^abs(r-3)
# φ₂₂₂(x,y,z)  = r*(r-1)*(r-2)*(x+y)^abs(r-3)
# φ₁₂₁(x,y,z)  = r*(r-1)*(r-2)*(x+y)^abs(r-3)
# φ₁₂₂(x,y,z)  = r*(r-1)*(r-2)*(x+y)^abs(r-3)

M₁₁(x, y, z) = -Dᵇ * (φ₁₁(x, y, z) + ν * φ₂₂(x, y, z))
M₁₂(x, y, z) = -Dᵇ * (1 - ν) * 0.5 * (φ₁₂(x, y, z) + φ₂₁(x, y, z))
M₂₂(x, y, z) = -Dᵇ * (ν * φ₁₁(x, y, z) + φ₂₂(x, y, z))
M₁₁₁(x, y, z) = -Dᵇ * (φ₁₁₁(x, y, z) + ν * φ₂₂₁(x, y, z))
M₁₂₂(x, y, z) = -Dᵇ * (1 - ν) * φ₁₂₂(x, y, z)
M₁₂₁(x, y, z) = -Dᵇ * (1 - ν) * φ₁₂₁(x, y, z)
M₂₂₂(x, y, z) = -Dᵇ * (ν * φ₁₁₂(x, y, z) + φ₂₂₂(x, y, z))

Q₁(x, y, z) = Dˢ * (w₁(x, y, z) - φ₁(x, y, z))
Q₂(x, y, z) = Dˢ * (w₂(x, y, z) - φ₂(x, y, z))
Q₁₁(x, y, z) = Dˢ * (w₁₁(x, y, z) - φ₁₁(x, y, z))
Q₂₂(x, y, z) = Dˢ * (w₂₂(x, y, z) - φ₂₂(x, y, z))
q(x, y, z) = -Q₁₁(x, y, z) - Q₂₂(x, y, z)
m₁(x, y, z) = M₁₁₁(x, y, z) + M₁₂₂(x, y, z) - Q₁(x, y, z)
m₂(x, y, z) = M₁₂₁(x, y, z) + M₂₂₂(x, y, z) - Q₂(x, y, z)
# ──────────────────────────────────────────────────────────

const to = TimerOutput()

gmsh.initialize()
# ──────────────────────────────────────────────────────────
type_w = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_φ = :(ReproducingKernel{:Linear2D,:□,:CubicSpline})
type_Q = :tri3
type_M = :(PiecewisePolynomial{:Linear2D})
ndiv_φ = 20
# ndiv_w = 33
# ndiv_q = 28
sʷ = 1.5
sᵠ = 1.5
for ndiv_w = 10:33
 XLSX.openxlsx("xls/square_eigen_$(ndiv_φ)_tri3_$(ndiv_w)_un.xlsx", mode="w") do xf
for ndiv_q = 10:33
row = ndiv_q
# ─── Deflection W ─────────────────────────────────────────
@timeit to "open msh file" gmsh.open("msh/patchtest_un_tri3_$ndiv_w.msh")
@timeit to "get nodes" nodes_w = get𝑿ᵢ()
@timeit to "get entities" entities_w = getPhysicalGroups()
xʷ = nodes_w.x
yʷ = nodes_w.y
zʷ = nodes_w.z
sp_w = RegularGrid(xʷ,yʷ,zʷ,n = 3,γ = 5)
nʷ = length(nodes_w)
# s = 1/ndiv_w
# s₁ = sʷ * s * ones(nʷ)
# s₂ = sʷ * s * ones(nʷ)
# s₃ = sʷ * s * ones(nʷ)
# push!(nodes_w,:s₁=>s₁,:s₂=>s₂,:s₃=>s₃)
elements_support = getElements(nodes_w, entities_w["Ω"], 1)
s_w, var_A = cal_area_support(elements_support)
s = sʷ*s_w*ones(nʷ)
push!(nodes_w, :s₁=>s, :s₂=>s, :s₃=>s)
# ─── Rotation Φ ───────────────────────────────────────────
@timeit to "open msh file" gmsh.open("msh/patchtest_un_tri3_$ndiv_φ.msh")
@timeit to "get nodes" nodes_φ = get𝑿ᵢ()
@timeit to "get entities" entities_φ = getPhysicalGroups()
xᵠ = nodes_φ.x
yᵠ = nodes_φ.y
zᵠ = nodes_φ.z
sp_φ = RegularGrid(xᵠ,yᵠ,zᵠ,n = 3,γ = 5)
nᵠ = length(nodes_φ)
# s = 1/ndiv_φ
# s₁ = sᵠ * s * ones(nᵠ)
# s₂ = sᵠ * s * ones(nᵠ)
# s₃ = sᵠ * s * ones(nᵠ)
# push!(nodes_φ,:s₁=>s₁,:s₂=>s₂,:s₃=>s₃)
elements_support = getElements(nodes_φ, entities_φ["Ω"], 1)
s_φ, var_A = cal_area_support(elements_support)
s = sᵠ*s_φ*ones(nᵠ)
push!(nodes_φ, :s₁=>s, :s₂=>s, :s₃=>s)
# ─── Shear ────────────────────────────────────────────────
@timeit to "open msh file" gmsh.open("msh/patchtest_un_tri3_$ndiv_q.msh")
@timeit to "get nodes" nodes = get𝑿ᵢ()
@timeit to "get entities" entities = getPhysicalGroups()

nˢ = length(nodes)
kᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)
kʷʷ = zeros(nʷ, nʷ)
kˢˢ = zeros(2 * nˢ, 2 * nˢ)
kˢᵠ = zeros(2 * nˢ, 2 * nᵠ)
kˢʷ = zeros(2 * nˢ, nʷ)
kᵠʷ = zeros(2 * nᵠ, nʷ)
fᵠ = zeros(2 * nᵠ)
fʷ = zeros(nʷ)
fˢ = zeros(2 * nˢ)
kᵅᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)
kᵅʷʷ = zeros(nʷ, nʷ)
fᵅᵠ = zeros(2 * nᵠ)
fᵅʷ = zeros(nʷ)

@timeit to "calculate ∫QQdΩ ∫∇QwdΩ" begin
    @timeit to "get elements" elements_q = getElements(nodes, entities["Ω"], integrationOrder)
    prescribe!(elements_q, :E=>E, :ν=>ν, :h=>h)
    @timeit to "calculate shape functions" set∇𝝭!(elements_q)

    @timeit to "get elements" elements_w = getElements(nodes_w, entities["Ω"], eval(type_w), integrationOrder, sp_w)
    prescribe!(elements_w, :E=>E, :ν=>ν, :h=>h, :q=>q)
    @timeit to "calculate shape functions" set𝝭!(elements_w)

    @timeit to "get elements" elements_w_Γ = getElements(nodes_w, entities["Γ"], eval(type_w), integrationOrder, sp_w, normal=true)
    @timeit to "calculate shape functions" set𝝭!(elements_w_Γ)

    @timeit to "get elements" elements_q_Γ = getElements(nodes, entities["Γ"], integrationOrder, normal=true)
    @timeit to "calculate shape functions" set𝝭!(elements_q_Γ)
 
    𝑎ˢˢ = ∫QQdΩ=>elements_q
    𝑎ˢʷ = [
        ∫∇QwdΩ=>(elements_q,elements_w),
        ∫QwdΓ=>(elements_q_Γ,elements_w_Γ),
    ]
    𝑓ʷ = ∫wqdΩ=>elements_w
    @timeit to "assemble" 𝑎ˢˢ(kˢˢ)
    @timeit to "assemble" 𝑎ˢʷ(kˢʷ)
    @timeit to "assemble" 𝑓ʷ(fʷ)
end

nₑ = length(elements_q)
nᵐ = nₑ*ApproxOperator.get𝑛𝑝(eval(type_M)(𝑿ᵢ[],𝑿ₛ[]))
kᵐᵐ = zeros(3*nᵐ,3*nᵐ)
kᵐᵠ = zeros(3*nᵐ,2*nᵠ)
kᵐʷ = zeros(3*nᵐ,nʷ)
kˢᵐ = zeros(2*nˢ,3*nᵐ)
fᵐ = zeros(3*nᵐ)

@timeit to "calculate ∫MMdΩ ∫MφdΩ" begin
    @timeit to "get elements" elements_m = getPiecewiseElements(entities["Ω"], eval(type_M), integrationOrder)
    prescribe!(elements_m, :E => E, :ν => ν, :h => h, :m₁ => m₁, :m₂ => m₂)
    @timeit to "calculate shape functions" set∇𝝭!(elements_m)

    @timeit to "get elements" elements_φ = getElements(nodes_φ, entities["Ω"], eval(type_φ), integrationOrder, sp_φ)
    prescribe!(elements_φ, :E => E, :ν => ν, :h => h, :m₁ => m₁, :m₂ => m₂)
    @timeit to "calculate shape functions" set𝝭!(elements_φ)

    @timeit to "get elements" elements_φ_Γ = getElements(nodes_φ, entities["Γ"], eval(type_φ), integrationOrder, sp_φ, normal=true)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_Γ)

    @timeit to "get elements" elements_m_Γ = getPiecewiseBoundaryElements(entities["Γ"], entities["Ω"], eval(type_M), integrationOrder)
    @timeit to "calculate shape functions" set𝝭!(elements_m_Γ)

    𝑎ᵐᵐ = ∫MMdΩ=>elements_m
    𝑎ᵐᵠ = [
        ∫∇MφdΩ=>(elements_m,elements_φ),
        ∫MφdΓ=>(elements_m_Γ,elements_φ_Γ),
    ]
    𝑎ˢᵠ = ∫QφdΩ=>(elements_q,elements_φ)
    𝑓ᵠ = ∫φmdΩ=>elements_φ
    @timeit to "assemble" 𝑎ᵐᵐ(kᵐᵐ)
    @timeit to "assemble" 𝑎ᵐᵠ(kᵐᵠ)
    @timeit to "assemble" 𝑎ˢᵠ(kˢᵠ)
    @timeit to "assemble" 𝑓ᵠ(fᵠ)
end

αʷ = 0e1;αᵠ = 0e1;
@timeit to "calculate ∫QwdΓ" begin
    @timeit to "get elements" elements_q_1 = getElements(nodes, entities["Γ¹"], integrationOrder, normal=true)
    @timeit to "get elements" elements_q_2 = getElements(nodes, entities["Γ²"], integrationOrder, normal=true)
    @timeit to "get elements" elements_q_3 = getElements(nodes, entities["Γ³"], integrationOrder, normal=true)
    @timeit to "get elements" elements_q_4 = getElements(nodes, entities["Γ⁴"], integrationOrder, normal=true)
    @timeit to "get elements" elements_w_1 = getElements(nodes_w, entities["Γ¹"], eval(type_w), integrationOrder, sp_w, normal=true)
    @timeit to "get elements" elements_w_2 = getElements(nodes_w, entities["Γ²"], eval(type_w), integrationOrder, sp_w, normal=true)
    @timeit to "get elements" elements_w_3 = getElements(nodes_w, entities["Γ³"], eval(type_w), integrationOrder, sp_w, normal=true)
    @timeit to "get elements" elements_w_4 = getElements(nodes_w, entities["Γ⁴"], eval(type_w), integrationOrder, sp_w, normal=true)
    prescribe!(elements_w_1, :α=>αʷ, :g=>w)
    prescribe!(elements_w_2, :α=>αʷ, :g=>w)
    prescribe!(elements_w_3, :α=>αʷ, :g=>w)
    prescribe!(elements_w_4, :α=>αʷ, :g=>w)
    @timeit to "calculate shape functions" set𝝭!(elements_q_1)
    @timeit to "calculate shape functions" set𝝭!(elements_q_2)
    @timeit to "calculate shape functions" set𝝭!(elements_q_3)
    @timeit to "calculate shape functions" set𝝭!(elements_q_4)
    @timeit to "calculate shape functions" set𝝭!(elements_w_1)
    @timeit to "calculate shape functions" set𝝭!(elements_w_2)
    @timeit to "calculate shape functions" set𝝭!(elements_w_3)
    @timeit to "calculate shape functions" set𝝭!(elements_w_4)
        𝑎 = ∫QwdΓ => (elements_q_1 ∪ elements_q_2 ∪ elements_q_3 ∪ elements_q_4, elements_w_1 ∪ elements_w_2 ∪ elements_w_3 ∪ elements_w_4)
    𝑎ʷ = ∫αwwdΓ => elements_w_1 ∪ elements_w_2 ∪ elements_w_3 ∪ elements_w_4
    @timeit to "assemble" 𝑎(kˢʷ,fˢ)
    @timeit to "assemble" 𝑎ʷ(kᵅʷʷ, fᵅʷ)
end

@timeit to "calculate ∫MφdΓ" begin
    @timeit to "get elements" elements_m_1 = getElements(entities["Γ¹"], entities["Γ"], elements_m_Γ)
    @timeit to "get elements" elements_m_2 = getElements(entities["Γ²"], entities["Γ"], elements_m_Γ)
    @timeit to "get elements" elements_m_3 = getElements(entities["Γ³"], entities["Γ"], elements_m_Γ)
    @timeit to "get elements" elements_m_4 = getElements(entities["Γ⁴"], entities["Γ"], elements_m_Γ)
    @timeit to "get elements" elements_φ_1 = getElements(nodes_φ, entities["Γ¹"], eval(type_φ), integrationOrder, sp_φ, normal=true)
    @timeit to "get elements" elements_φ_2 = getElements(nodes_φ, entities["Γ²"], eval(type_φ), integrationOrder, sp_φ, normal=true)
    @timeit to "get elements" elements_φ_3 = getElements(nodes_φ, entities["Γ³"], eval(type_φ), integrationOrder, sp_φ, normal=true)
    @timeit to "get elements" elements_φ_4 = getElements(nodes_φ, entities["Γ⁴"], eval(type_φ), integrationOrder, sp_φ, normal=true)
    prescribe!(elements_φ_1, :α=>αᵠ, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_φ_2, :α=>αᵠ, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_φ_3, :α=>αᵠ, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_φ_4, :α=>αᵠ, :g₁=>φ₁, :g₂=>φ₂, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_1)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_2)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_3)
    @timeit to "calculate shape functions" set𝝭!(elements_φ_4)
    𝑎 = ∫MφdΓ => (elements_m_1 ∪ elements_m_2 ∪ elements_m_3 ∪ elements_m_4, elements_φ_1 ∪ elements_φ_2 ∪ elements_φ_3 ∪ elements_φ_4)
    𝑎ᵅ = ∫αφφdΓ => elements_φ_1 ∪ elements_φ_2 ∪ elements_φ_3 ∪ elements_φ_4
    @timeit to "assemble" 𝑎(kᵐᵠ,fᵐ)
    @timeit to "assemble" 𝑎ᵅ(kᵅᵠᵠ, fᵅᵠ)
end

# @timeit to "solve" d = [kᵠᵠ+kᵅᵠᵠ kᵠʷ kˢᵠ' kᵐᵠ';kᵠʷ' kʷʷ+kᵅʷʷ kˢʷ' kᵐʷ';kˢᵠ kˢʷ kˢˢ kˢᵐ;kᵐᵠ kᵐʷ kˢᵐ' kᵐᵐ]\[fᵠ+fᵅᵠ;fʷ+fᵅʷ;fˢ;fᵐ]
@timeit to "calculate ∫wwdΩ" begin
    @timeit to "get elements" elements = getElements(nodes_w, entities["Ω"], eval(type_w), integrationOrder, sp_w)
    prescribe!(elements, :E=>E, :ν=>ν, :h=>h)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    𝑎ʷʷ = ∫wwdΩ=>elements
    @timeit to "assemble" 𝑎ʷʷ(kʷʷ)
end



println("h = $h, Dˢ = $Dˢ, Dᵇ = $Dᵇ, nᵠ = $nᵠ, nʷ = $nʷ, nˢ = $nˢ")
print("nˢ≤ᵠ:         ")
n_diff = nᵠ-nˢ
n_diff≥0.0 ? println("✓:$n_diff") : println("×:$n_diff")
print("nʷ≤⌊[nˢ]⌋-1:  ")
n = floor(0.5*((1+8*nˢ)^0.5-3))
n_diff = 0.5*n*(n+1)-nʷ
n_diff≥0.0 ? println("✓:$n_diff") : println("×:$n_diff")

# ─── Eigen Test For βʷ ────────────────────────────────────
βʷ² = eigvals(-kˢʷ'*(kˢˢ\kˢʷ))
# βʷ² = eigvals(kᵠʷ*(kʷʷ\kᵠʷ'), kᵠᵠ)
# βʷ² = eigvals(-kˢʷ*(kʷʷ\kˢʷ'),kˢˢ)
βʷ² = real.(βʷ²)
βʷ²⁺ = βʷ²[βʷ² .≥ 1e6*eps()]
βʷ⁺ = βʷ²⁺.^0.5
nʷ⁺ = length(βʷ⁺)
# println(βʷ⁺)
βʷ⁺ = min(βʷ⁺...)
println("βʷ⁺ = $βʷ⁺, nʷ⁺ = $nʷ⁺")

# ─── Eigen Test For βᵞ ────────────────────────────────────
k̃ᵠᵠ = - kᵐᵠ'*(kᵐᵐ\kᵐᵠ)
k̃ʷʷ = - kˢʷ'*(kˢˢ\kˢʷ)
k̃ᵠʷ = - kˢᵠ'*(kˢˢ\kˢʷ)
βᵞ² = eigvals(-kˢᵠ'*(kˢˢ\kˢᵠ)-k̃ᵠʷ*(k̃ʷʷ\k̃ᵠʷ'),k̃ᵠᵠ)
# βᵞ² = eigvals(k̃ᵠᵠ)
βᵞ² = real.(βᵞ²)
βᵞ²⁺ = βᵞ²[βᵞ² .≥ 1e6*eps()]
βᵞ⁺ = βᵞ²⁺.^0.5
nᵞ⁺ = length(βᵞ⁺)
βᵞ⁺ = min(βᵞ⁺...)
println("βᵞ⁺ = $βᵞ⁺, nᵞ⁺ = $nᵞ⁺")
# ──────────────────────────────────────────────────────────
println("L₂ error of βʷ⁺: ", log10(βʷ⁺))
println("L₂ error of βᵞ⁺: ", log10(βᵞ⁺))
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
     sheet["I1"] = "βʷ⁺"
     sheet["J1"] = "βᵞ⁺"
     sheet["K1"] = "L₂βʷ⁺"
     sheet["L1"] = "L₂βᵞ⁺"

     sheet["A$row"] = "$type_w"
     sheet["B$row"] = nʷ
     sheet["C$row"] = "$type_φ"
     sheet["D$row"] = nᵠ
     sheet["E$row"] = "$type_Q"
     sheet["F$row"] = nˢ
     sheet["G$row"] = "$type_M"
     sheet["H$row"] = nᵐ
     sheet["I$row"] = "$βʷ⁺"
     sheet["J$row"] = "$βᵞ⁺"
     sheet["K$row"] = log10(βʷ⁺)
     sheet["L$row"] = log10(βᵞ⁺)

end
end
end

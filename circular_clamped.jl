using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, L₂, L₂φ, ∫wVdΓ, ∫φMdΓ

using LinearAlgebra
using TimerOutputs
import Gmsh: gmsh

# -----------------------------
# case switch (reference solution from 圆板精確解.md)
# -----------------------------

    E = 10.92
    ν = 0.3
    h = 1.0
    R = 5.0
    fz = 1.0
    α = 1e8 * E 

    Dᵇ = E*h^3/12/(1-ν^2)
    Dˢ = 5/6*E*h/(2*(1+ν))
    
 

# -----------------------------
# r = 1

# w(x, y, z) = fz * R^4 / (64 * Dᵇ) * (1 - ξ(x, y)^2) * ((1 - ξ(x, y)^2) + (8 * (h / R)^2) / (3 * k * (1 - ν)))

# # w(x, y, z) = (x + y)^r
# w₁(x, y, z) = r * (x + y)^abs(r - 1)
# w₂(x, y, z) = r * (x + y)^abs(r - 1)
# w₁₁(x, y, z) = r * (r - 1) * (x + y)^abs(r - 2)
# w₂₂(x, y, z) = r * (r - 1) * (x + y)^abs(r - 2)
# φ₁(x, y, z) = r * (x + y)^abs(r - 1)
# φ₂(x, y, z) = r * (x + y)^abs(r - 1)
# φ₁₁(x, y, z) = r * (r - 1) * (x + y)^abs(r - 2)
# φ₁₂(x, y, z) = r * (r - 1) * (x + y)^abs(r - 2)
# φ₂₁(x, y, z) = r * (r - 1) * (x + y)^abs(r - 2)
# φ₂₂(x, y, z) = r * (r - 1) * (x + y)^abs(r - 2)
# φ₁₁₁(x, y, z) = r * (r - 1) * (r - 2) * (x + y)^abs(r - 3)
# φ₁₁₂(x, y, z) = r * (r - 1) * (r - 2) * (x + y)^abs(r - 3)
# φ₂₂₁(x, y, z) = r * (r - 1) * (r - 2) * (x + y)^abs(r - 3)
# φ₂₂₂(x, y, z) = r * (r - 1) * (r - 2) * (x + y)^abs(r - 3)
# φ₁₂₁(x, y, z) = r * (r - 1) * (r - 2) * (x + y)^abs(r - 3)
# φ₁₂₂(x, y, z) = r * (r - 1) * (r - 2) * (x + y)^abs(r - 3)

# # M₁₁(x, y, z) = -Dᵇ * (φ₁₁(x, y, z) + ν * φ₂₂(x, y, z))
# # M₁₂(x, y, z) = -Dᵇ * (1 - ν) * 0.5 * (φ₁₂(x, y, z) + φ₂₁(x, y, z))
# # M₂₂(x, y, z) = -Dᵇ * (ν * φ₁₁(x, y, z) + φ₂₂(x, y, z))
# M₁₁₁(x, y, z) = -Dᵇ * (φ₁₁₁(x, y, z) + ν * φ₂₂₁(x, y, z))
# M₁₂₂(x, y, z) = -Dᵇ * (1 - ν) * φ₁₂₂(x, y, z)
# M₁₂₁(x, y, z) = -Dᵇ * (1 - ν) * φ₁₂₁(x, y, z)
# M₂₂₂(x, y, z) = -Dᵇ * (ν * φ₁₁₂(x, y, z) + φ₂₂₂(x, y, z))

# Q₁(x, y, z) = Dˢ * (w₁(x, y, z) - φ₁(x, y, z))
# Q₂(x, y, z) = Dˢ * (w₂(x, y, z) - φ₂(x, y, z))
# Q₁₁(x, y, z) = Dˢ * (w₁₁(x, y, z) - φ₁₁(x, y, z))
# Q₂₂(x, y, z) = Dˢ * (w₂₂(x, y, z) - φ₂₂(x, y, z))
# q(x, y, z) = -Q₁₁(x, y, z) - Q₂₂(x, y, z)
# m₁(x, y, z) = M₁₁₁(x, y, z) + M₁₂₂(x, y, z) - Q₁(x, y, z)
# m₂(x, y, z) = M₁₂₁(x, y, z) + M₂₂₂(x, y, z) - Q₂(x, y, z)
# # ---------------------------------------------------------------------------------------------------
ξ(x, y) = sqrt(x^2 + y^2) / R
Mr(x, y, z) = fz * R^2 / 16 * (1 + ν) * (1 - ((3 + ν) / (1 + ν)) * ξ(x, y)^2)
Mθ(x, y, z) = fz * R^2 / 16 * (1 + ν) * (1 - ((1 + 3ν) / (1 + ν)) * ξ(x, y)^2)
# # Pi_int() = (fz^2 * R^6 * ν) / (384 * Dᵇ) * (1 + (4 * (h / R)^2) / (k * (1 - ν)))
Tr(x, y, z) = -fz * axis_r(x, y) / 2
# ---------------------------------------------------------------------------------------------------


# 自然
# （示意图边界验收时用零边界；解析解对比时也依然是 clamped 零边界）
# w0(x, y, z) = 0.0
# φ10(x, y, z) = 0.0
# φ20(x, y, z) = 0.0

# V(x, y, z) = 0.0
# M₁(x, y, z) = 0.0
# M₂(x, y, z) = 0.0



const to = TimerOutput()

gmsh.initialize()
@timeit to "open msh file" gmsh.open("./msh/circular.msh")
# @timeit to "open msh file" gmsh.open("./msh/circular_tri3_9.msh")
# @timeit to "open msh file" gmsh.open("./msh/circular_tri3_15.msh")
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()

nʷ = length(nodes)
nᵠ = length(nodes)

kʷʷ = zeros(nʷ, nʷ)
kᵠᵠ = zeros(2 * nᵠ, 2 * nᵠ)
kᵠʷ = zeros(2 * nᵠ, nʷ)
fʷ = zeros(nʷ)
fᵠ = zeros(2 * nᵠ)

@timeit to "assemble domain" begin
    elements = getElements(nodes, entities["Ω"])
    prescribe!(elements, :E => E, :ν => ν, :h => h, :q => fz)
    set∇𝝭!(elements)

    𝑎ʷʷ = ∫wwdΩ => elements
    𝑎ᵠʷ = ∫φwdΩ => elements
    𝑎ᵠᵠ = [
        ∫φφdΩ => elements,
        ∫κκdΩ => elements,
    ]
    𝑓ʷ = ∫wqdΩ => elements
    # 𝑓ᵠ = ∫φmdΩ => elements
    @timeit to "assemble" 𝑎ʷʷ(kʷʷ)
    @timeit to "assemble" 𝑎ᵠʷ(kᵠʷ)
    @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)
    @timeit to "assemble" 𝑓ʷ(fʷ)
    # @timeit to "assemble" 𝑓ᵠ(fᵠ)
end

# println("after domain: ‖fʷ‖₂=", norm(fʷ), "  ‖fᵠ‖₂=", norm(fᵠ))

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
        :n₂₂ => 1.0
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
        :g => w,
        :g₁ => φ₁,
        :g₂ => φ₂,
        :n₁₁ => 1.0,
        :n₁₂ => 0.0,
        :n₂₂ => 0.0,
        # :V => (x,y,z,n₁,n₂)->Q₁(x,y,z)*n₁ + Q₂(x,y,z)*n₂,
        :M₁ => 0.0,
        :M₂ => (x,y,z,n₁,n₂)->Mθ(x, y, z)*n₁ + Mr(x, y, z)*n₂,


    )
    𝑎ʷ = ∫αwwdΓ => elements
    𝑎ᵠ = ∫αφφdΓ => elements
    # 𝑓ʷ = ∫wVdΓ => elements
    𝑓ᵠ = ∫φMdΓ => elements
    # @timeit to "assemble" 𝑎ʷ(kʷʷ,fʷ)
    @timeit to "assemble" 𝑎ᵠ(kᵠᵠ,fᵠ)
    @timeit to "assemble" 𝑓ʷ(fʷ)
    @timeit to "assemble" 𝑓ᵠ(fᵠ)
end

@timeit to "Γᵇ" begin
    elements = getElements(nodes, entities["Γᵇ"],normal=true)
    set𝝭!(elements)

    prescribe!(elements, 
        :α => α,
        :g => w,
        :g₁ => φ₁,
        :g₂ => φ₂,
        :n₁₁ => 0.0,
        :n₁₂ => 0.0,
        :n₂₂ => 1.0,
        # :V => (x,y,z,n₁,n₂)->Q₁(x,y,z)*n₁ + Q₂(x,y,z)*n₂,
        :M₁ => (x,y,z,n₁,n₂)->Mr(x, y, z)*n₁ + Mθ(x, y, z)*n₂,
        :M₂ => 0.0,


    )
    𝑎ʷ = ∫αwwdΓ => elements
    𝑎ᵠ = ∫αφφdΓ => elements
    𝑓ʷ = ∫wVdΓ => elements
    𝑓ᵠ = ∫φMdΓ => elements
    # @timeit to "assemble" 𝑎ʷ(kʷʷ,fʷ)
    @timeit to "assemble" 𝑎ᵠ(kᵠᵠ,fᵠ)
    # @timeit to "assemble" 𝑓ʷ(fʷ)
    @timeit to "assemble" 𝑓ᵠ(fᵠ)
end

# println("after bnd: ‖fʷ‖₂=", norm(fʷ), "  ‖fᵠ‖₂=", norm(fᵠ))

@timeit to "solve" d = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ] \ [fᵠ; fʷ]

# println("check d norms: ‖d‖₂=", norm(d), "  ‖f‖₂=", norm([fᵠ; fʷ]))
# println("check rhs norms: ‖fʷ‖₂=", norm(fʷ), "  ‖fᵠ‖₂=", norm(fᵠ))

push!(nodes, :d => d[2*nᵠ+1:end], :d₁ => d[1:2:2*nᵠ], :d₂ => d[2:2:2*nᵠ])

# println("check node fields:")
# println("  node1: d=", nodes[1].d, " d₁=", nodes[1].d₁, " d₂=", nodes[1].d₂)
# println("  node2: d=", nodes[2].d, " d₁=", nodes[2].d₁, " d₂=", nodes[2].d₂)


@timeit to "calculate error" begin
    elements_err = getElements(nodes, entities["Ω"], 10)
    prescribe!(elements_err, :E => E, :ν => ν, :h => h, :u => w, :φ₁ => φ₁, :φ₂ => φ₂)
    set𝝭!(elements_err)
    global L₂_w = L₂(elements_err)
    global L₂_φ = L₂φ(elements_err)
end


gmsh.finalize()

println(to)

println("α penalty: ", α)
println("L₂ error of w: ", L₂_w)
println("L₂ error of φ: ", L₂_φ)

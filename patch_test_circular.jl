using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, L₂, L₂φ, ∫wVdΓ, ∫φMdΓ

using TimerOutputs
import Gmsh: gmsh

# penalty（与 patch_test.jl 一致）

E = 1.0
ν = 0.3
h = 1e-1
Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))
α = 1e8 * E

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
@timeit to "open msh file" gmsh.open("./msh/circular_tri3_3.msh")
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
    prescribe!(elements, :E => E, :ν => ν, :h => h, :q => q, :m₁ => m₁, :m₂ => m₂)
    set∇𝝭!(elements)

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
        :V => (x,y,z,n₁,n₂)->Q₁(x,y,z)*n₁ + Q₂(x,y,z)*n₂,
        :M₁ => 0.0,
        :M₂ => (x,y,z,n₁,n₂)->M₁₂(x,y,z)*n₁ + M₂₂(x,y,z)*n₂,
    )
    𝑎ʷ = ∫αwwdΓ => elements
    𝑎ᵠ = ∫αφφdΓ => elements
    𝑓ʷ = ∫wVdΓ => elements
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
        :V => (x,y,z,n₁,n₂)->Q₁(x,y,z)*n₁ + Q₂(x,y,z)*n₂,
        :M₁ => (x,y,z,n₁,n₂)->M₁₁(x,y,z)*n₁ + M₁₂(x,y,z)*n₂,
        :M₂ => 0.0,
    )
    𝑎ʷ = ∫αwwdΓ => elements
    𝑎ᵠ = ∫αφφdΓ => elements
    𝑓ʷ = ∫wVdΓ => elements
    𝑓ᵠ = ∫φMdΓ => elements
    # @timeit to "assemble" 𝑎ʷ(kʷʷ,fʷ)
    @timeit to "assemble" 𝑎ᵠ(kᵠᵠ,fᵠ)
    @timeit to "assemble" 𝑓ʷ(fʷ)
    @timeit to "assemble" 𝑓ᵠ(fᵠ)
end

@timeit to "solve" d = [kᵠᵠ kᵠʷ; kᵠʷ' kʷʷ] \ [fᵠ; fʷ]
push!(nodes, :d => d[2*nᵠ+1:end], :d₁ => d[1:2:2*nᵠ], :d₂ => d[2:2:2*nᵠ])

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
using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫wVdΓ, ∫φMdΓ, ∫αwwdΓ, ∫αφφdΓ, L₂, L₂φ

using TimerOutputs 
import Gmsh: gmsh

E = 1.0
ν = 0.3
h = 10^-1
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

gmsh.initialize()
@timeit to "open msh file" gmsh.open("msh/patchtest.msh")
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
end

@timeit to "calculate ∫αwwdΓ ∫αφφdΓ" begin
    @timeit to "get elements" elements_1 = getElements(nodes, entities["Γ¹"])
    @timeit to "get elements" elements_2 = getElements(nodes, entities["Γ²"])
    @timeit to "get elements" elements_3 = getElements(nodes, entities["Γ³"])
    @timeit to "get elements" elements_4 = getElements(nodes, entities["Γ⁴"])
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
    @timeit to "assemble" 𝑎ʷ(kʷʷ,fʷ)
    @timeit to "assemble" 𝑎ᵠ(kᵠᵠ,fᵠ)
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




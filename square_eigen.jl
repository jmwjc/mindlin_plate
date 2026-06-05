using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫wVdΓ, ∫φMdΓ, ∫αwwdΓ, ∫αφφdΓ

using TimerOutputs, LinearAlgebra
import Gmsh: gmsh

E = 1e6
ν = 0.3
h = 1e-3
Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))
integrationOrder = 2

const to = TimerOutput()

gmsh.initialize()
@timeit to "open msh file" gmsh.open("msh/patchtest_quad4_32.msh")
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
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"],integrationOrder)
    @timeit to "get elements" elements_reduced = getElements(nodes, entities["Ω"],1)
    prescribe!(elements, :E=>E, :ν=>ν, :h=>h)
    prescribe!(elements_reduced, :E=>E, :ν=>ν, :h=>h)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    @timeit to "calculate shape functions" set∇𝝭!(elements_reduced)
    𝑎ʷʷ = ∫wwdΩ=>elements_reduced
    𝑎ᵠʷ = ∫φwdΩ=>elements_reduced
    𝑎ᵠᵠ = [
        ∫φφdΩ=>elements_reduced,
        ∫κκdΩ=>elements,
    ]
    @timeit to "assemble" 𝑎ʷʷ(kʷʷ)
    @timeit to "assemble" 𝑎ᵠʷ(kᵠʷ)
    @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)
end

@timeit to "calculate ∫αwwdΓ ∫αφφdΓ" begin
    @timeit to "get elements" elements_1 = getElements(nodes, entities["Γ¹"],integrationOrder)
    @timeit to "get elements" elements_2 = getElements(nodes, entities["Γ²"],integrationOrder)
    @timeit to "get elements" elements_3 = getElements(nodes, entities["Γ³"],integrationOrder)
    @timeit to "get elements" elements_4 = getElements(nodes, entities["Γ⁴"],integrationOrder)
    prescribe!(elements_1, :α=>1e8*E, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_2, :α=>1e8*E, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_3, :α=>1e8*E, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_4, :α=>1e8*E, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    @timeit to "calculate shape functions" set𝝭!(elements_4)
    𝑎ʷ = ∫αwwdΓ=>elements_1∪elements_2∪elements_3∪elements_4
    𝑎ᵠ = ∫αφφdΓ=>elements_1∪elements_2∪elements_3∪elements_4
    @timeit to "assemble" 𝑎ʷ(kʷʷ)
    @timeit to "assemble" 𝑎ᵠ(kᵠᵠ)
end


println("h = $h, Dˢ = $Dˢ, Dᵇ = $Dᵇ, nᵠ = $nᵠ, nʷ = $nʷ")
β² = eigvals(kᵠʷ*(kʷʷ\kᵠʷ'), kᵠᵠ)
β² = real.(β²)
β²⁺ = β²[β² .≥ 1e4*eps()]
β⁺ = β²⁺.^0.5
n⁺ = length(β⁺)
println("β⁺ = $β⁺, n⁺ = $n⁺")
# β² = eigvals(kᵠᵠ, -kᵠʷ*(kʷʷ\kᵠʷ'))
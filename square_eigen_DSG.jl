using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ, ∫φφdΩ_DSG, ∫φwdΩ_DSG, ∫wqdΩ, ∫φmdΩ, ∫wVdΓ, ∫φMdΓ, ∫αwwdΓ, ∫αφφdΓ, L₂w, L₂φ

using TimerOutputs, LinearAlgebra 
import Gmsh: gmsh

E = 10.92e6
ν = 0.3
h = 1e-3
Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))

const to = TimerOutput()

ndiv = 4
integrationOrder = 2

gmsh.initialize()
@timeit to "open msh file" gmsh.open("msh/patchtest_tri3_$ndiv.msh")
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get nodes" nodes = get𝑿ᵢ()

nʷ = length(nodes)
nᵠ = length(nodes)
kʷʷ = zeros(nʷ,nʷ)
kᵠᵠᵇ = zeros(2*nᵠ,2*nᵠ)
kᵠᵠˢ = zeros(2*nᵠ,2*nᵠ)
kᵠʷ = zeros(2*nᵠ,nʷ)
fʷ = zeros(nʷ)
fᵠ = zeros(2*nᵠ)

@timeit to "calculate ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫wφdΩ" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"],integrationOrder)
    prescribe!(elements, :E=>E, :ν=>ν, :h=>h)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    𝑎ʷʷ = ∫wwdΩ=>elements
    𝑎ᵠʷ = ∫φwdΩ_DSG=>elements
    𝑎ᵠᵠˢ = ∫φφdΩ_DSG=>elements
    𝑎ᵠᵠᵇ = ∫κκdΩ=>elements
    @timeit to "assemble" 𝑎ʷʷ(kʷʷ)
    @timeit to "assemble" 𝑎ᵠʷ(kᵠʷ)
    @timeit to "assemble" 𝑎ᵠᵠˢ(kᵠᵠˢ)
    @timeit to "assemble" 𝑎ᵠᵠᵇ(kᵠᵠᵇ)
end

α = 1e8
@timeit to "calculate ∫αwwdΓ ∫αφφdΓ" begin
    @timeit to "get elements" elements_1 = getElements(nodes, entities["Γ¹"],integrationOrder)
    @timeit to "get elements" elements_2 = getElements(nodes, entities["Γ²"],integrationOrder)
    @timeit to "get elements" elements_3 = getElements(nodes, entities["Γ³"],integrationOrder)
    @timeit to "get elements" elements_4 = getElements(nodes, entities["Γ⁴"],integrationOrder)
    prescribe!(elements_1, :α=>α*E, :g=>0.0, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_2, :α=>α*E, :g=>0.0, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_3, :α=>α*E, :g=>0.0, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    prescribe!(elements_4, :α=>α*E, :g=>0.0, :g₁=>0.0, :g₂=>0.0, :n₁₁=>1.0, :n₁₂=>0.0, :n₂₂=>1.0)
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    @timeit to "calculate shape functions" set𝝭!(elements_4)
    𝑎ʷ = ∫αwwdΓ=>elements_1∪elements_2∪elements_3∪elements_4
    𝑎ᵠ = ∫αφφdΓ=>elements_1∪elements_2∪elements_3∪elements_4
    @timeit to "assemble" 𝑎ʷ(kʷʷ)
    @timeit to "assemble" 𝑎ᵠ(kᵠᵠᵇ)
end

gmsh.finalize()

# β² = eigvals(kᵠᵠᵇ)
β² = eigvals(kᵠᵠˢ-kᵠʷ*(kʷʷ\kᵠʷ'), kᵠᵠᵇ)
# β² = eigvals([kᵠᵠˢ kᵠʷ;kᵠʷ' kʷʷ], [kᵠᵠᵇ zeros(2nᵠ,nʷ);zeros(nʷ,2nᵠ) zeros(nʷ,nʷ)])

β² = real.(β²)
β²⁺ = β²[β² .≥ 1e2*eps()]
β⁺ = β²⁺.^0.5
n⁺ = length(β⁺)
β⁺ = min(β⁺...)
logβ⁺ = log10(β⁺)
println("β⁺ = $β⁺, n⁺ = $n⁺, logβ⁺ = $logβ⁺")
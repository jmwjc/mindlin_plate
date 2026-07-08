using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫∇w∇wdΩ, ∫φφdΩ, ∫φwdΩ, ∫φφdΩ_DSG, ∫φwdΩ_DSG, ∫wqdΩ, ∫φmdΩ, ∫wVdΓ, ∫φMdΓ, ∫αwwdΓ, ∫αφφdΓ, L₂w, L₂φ

using TimerOutputs, LinearAlgebra 
import Gmsh: gmsh

E = 10.92e6
ν = 0.3
h = 1e-3
Dᵇ = E/12/(1-ν^2)
Dˢ = 5/6*E/h^2/(2*(1+ν))

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
kᵠᵠ = zeros(2*nᵠ,2*nᵠ)
kᵠʷ = zeros(2*nᵠ,nʷ)
fʷ = zeros(nʷ)
fᵠ = zeros(2*nᵠ)

@timeit to "calculate ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫wφdΩ" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"],integrationOrder)
    prescribe!(elements, :E=>E, :ν=>ν, :h=>h)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    𝑎ʷʷ = ∫∇w∇wdΩ=>elements
    𝑎ᵠʷ = ∫φwdΩ_DSG=>elements
    𝑎ᵠᵠ = ∫φφdΩ_DSG=>elements
    @timeit to "assemble" 𝑎ʷʷ(kʷʷ)
    @timeit to "assemble" 𝑎ᵠʷ(kᵠʷ)
    @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)
end

gmsh.finalize()

β² = eigvals([kᵠᵠ kᵠʷ;kᵠʷ' kʷʷ]/Dˢ*(1/ndiv)^(-2))

β² = real.(β²)
β²⁺ = β²[β² .≥ 1e4*eps()]
β⁺ = β²⁺.^0.5
n⁺ = length(β⁺)
β⁺ = min(β⁺...)
logβ⁺ = log10(β⁺)
println("β⁺ = $β⁺, n⁺ = $n⁺, logβ⁺ = $logβ⁺")
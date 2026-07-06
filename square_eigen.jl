using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫∇w∇wdΩ, ∫φφdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫wVdΓ, ∫φMdΓ, ∫αwwdΓ, ∫αφφdΓ

using TimerOutputs, LinearAlgebra
import Gmsh: gmsh

E = 10.92e6
ν = 0.3
h = 1e-3
Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))
integrationOrder = 2
ndiv = 16
const to = TimerOutput()

gmsh.initialize()
# @timeit to "open msh file" gmsh.open("msh/patchtest_quad4_$ndiv.msh")
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
    @timeit to "get elements" elements_reduced = getElements(nodes, entities["Ω"],2)
    prescribe!(elements, :E=>E, :ν=>ν, :h=>h)
    prescribe!(elements_reduced, :E=>E, :ν=>ν, :h=>h)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    @timeit to "calculate shape functions" set∇𝝭!(elements_reduced)
    𝑎ʷʷ = ∫∇w∇wdΩ=>elements_reduced
    𝑎ᵠʷ = ∫φwdΩ=>elements_reduced
    𝑎ᵠᵠ = [
        ∫φφdΩ=>elements_reduced,
        ∫κκdΩ=>elements,
    ]
    @timeit to "assemble" 𝑎ʷʷ(kʷʷ)
    @timeit to "assemble" 𝑎ᵠʷ(kᵠʷ)
    @timeit to "assemble" 𝑎ᵠᵠ(kᵠᵠ)
end


println("h = $h, Dˢ = $Dˢ, Dᵇ = $Dᵇ, nᵠ = $nᵠ, nʷ = $nʷ")
h = 1/ndiv
# β² = eigvals([kᵠᵠ kᵠʷ;kᵠʷ' kʷʷ])
β² = eigvals(kʷʷ/Dˢ/h^2)
β² = real.(β²)
β²⁺ = β²[β² .≥ 1e4*eps()]
β⁺ = β²⁺.^0.5
n⁺ = length(β⁺)
β⁺ = min(β⁺...)
println("β⁺ = $β⁺, n⁺ = $n⁺")
logβ⁺ = log10(β⁺)
println("$logβ⁺")
using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.Heat: ∫∫∇v∇udxdy, ∫vbdΩ, ∫vtdΓ, ∫vgdΓ, H₁

using TimerOutputs
import Gmsh: gmsh

ndiv = 32
integrationOrder = 2

r = 5
𝑢(x,y,z) = (x+y)^r
∂𝑢∂x(x,y,z) = r*(x+y)^abs(r-1)
∂𝑢∂y(x,y,z) = r*(x+y)^abs(r-1)
𝑏(x,y,z) = -2*r*(r-1)*(x+y)^abs(r-2)

const to = TimerOutput()

gmsh.initialize()

@timeit to "open msh file" gmsh.open("msh/patchtest_tri3"*"_$ndiv.msh")
@timeit to "get nodes" nodes = get𝑿ᵢ()
@timeit to "get entities" entities = getPhysicalGroups()
@timeit to "get elements" elements = getElements(nodes, entities["Ω"], integrationOrder)
nₚ = length(nodes)
k = zeros(nₚ,nₚ)
f = zeros(nₚ)

@timeit to "calculate ∫∫∇v∇udxdy" begin
    prescribe!(elements, :k=>1.0, :b=>𝑏)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    𝑎 = ∫∫∇v∇udxdy=>elements
    𝑓 = ∫vbdΩ=>elements

    @timeit to "assemble" 𝑎(k)
    @timeit to "assemble" 𝑓(f)
end

α = 1e8
@timeit to "calculate ∫vgdΓ" begin
    @timeit to "get elements" elements_1 = getElements(nodes, entities["Γ¹"], integrationOrder)
    @timeit to "get elements" elements_2 = getElements(nodes, entities["Γ²"], integrationOrder)
    @timeit to "get elements" elements_3 = getElements(nodes, entities["Γ³"], integrationOrder)
    @timeit to "get elements" elements_4 = getElements(nodes, entities["Γ⁴"], integrationOrder)
    @timeit to "calculate shape functions" set𝝭!(elements_1)
    @timeit to "calculate shape functions" set𝝭!(elements_2)
    @timeit to "calculate shape functions" set𝝭!(elements_3)
    @timeit to "calculate shape functions" set𝝭!(elements_4)
    prescribe!(elements_1, :g=>𝑢, :α=>α)
    prescribe!(elements_2, :g=>𝑢, :α=>α)
    prescribe!(elements_3, :g=>𝑢, :α=>α)
    prescribe!(elements_4, :g=>𝑢, :α=>α)
    𝑓 = ∫vgdΓ=>elements_1∪elements_2∪elements_3∪elements_4

    @timeit to "assemble" 𝑓(k,f)
end

d = k\f

@timeit to "calculate error" begin
    @timeit to "get elements" elements = getElements(nodes, entities["Ω"], 10)
    push!(nodes,:d=>d)
    prescribe!(elements, :k=>1.0, :u=>𝑢, :∂u∂x=>∂𝑢∂x, :∂u∂y=>∂𝑢∂y, :∂u∂z=>0.0)
    @timeit to "calculate shape functions" set∇𝝭!(elements)
    H₁_error, L₂_error = H₁(elements) 
end

h = log10(1/ndiv)
H₁_error = log10(H₁_error)
L₂_error = log10(L₂_error)
# println("h=$h, L₂=$L₂_error, H₁=$H₁_error")
# println("$h, $L₂_error, $H₁_error")
println("$L₂_error, $H₁_error")
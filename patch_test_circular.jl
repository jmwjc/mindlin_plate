using ApproxOperator
import ApproxOperator.GmshImport: getPhysicalGroups, get𝑿ᵢ, getElements
import ApproxOperator.MindlinPlate: ∫κκdΩ, ∫wwdΩ, ∫φφdΩ, ∫φwdΩ, ∫wqdΩ, ∫φmdΩ, ∫αwwdΓ, ∫αφφdΓ, L₂, L₂φ, ∫wVdΓ, ∫φMdΓ

using TimerOutputs
import Gmsh: gmsh

E = 1.0
ν = 0.3
h = 1e-1
Dᵇ = E * h^3 / 12 / (1 - ν^2)
Dˢ = 5 / 6 * E * h / (2 * (1 + ν))

# penalty（与 patch_test.jl 一致）
α = 1e8 * E

# 制造解（阶段一严格线性 patch：w=ax+by+c，φ=∇w 常数；所有高阶导数为 0）
w(x, y, z) = 1.0 + x + y
w₁(x, y, z) = 1.0
w₂(x, y, z) = 1.0
w₁₁(x, y, z) = 0.0
w₂₂(x, y, z) = 0.0

# φ=∇w
φ₁(x, y, z) = w₁(x, y, z)
φ₂(x, y, z) = w₂(x, y, z)

# φ 的各阶导数（没有的填 0）
φ₁₁(x, y, z) = 0.0
φ₁₂(x, y, z) = 0.0
φ₂₁(x, y, z) = 0.0
φ₂₂(x, y, z) = 0.0
φ₁₁₁(x, y, z) = 0.0
φ₁₁₂(x, y, z) = 0.0
φ₂₂₁(x, y, z) = 0.0
φ₂₂₂(x, y, z) = 0.0
φ₁₂₁(x, y, z) = 0.0
φ₁₂₂(x, y, z) = 0.0

# 对严格线性 patch：M=0, Q=0 ⇒ q=m=0（仍显式写全，便于后续升级阶次）
M₁₁(x, y, z) = 0.0
M₁₂(x, y, z) = 0.0
M₂₂(x, y, z) = 0.0
M₁₁₁(x, y, z) = 0.0
M₁₂₂(x, y, z) = 0.0
M₁₂₁(x, y, z) = 0.0
M₂₂₂(x, y, z) = 0.0

Q₁(x, y, z) = 0.0
Q₂(x, y, z) = 0.0
Q₁₁(x, y, z) = 0.0
Q₂₂(x, y, z) = 0.0

q(x, y, z) = 0.0
m₁(x, y, z) = 0.0
m₂(x, y, z) = 0.0

# 边界项写全：阶段一先全部填 0（之后替换成对应精确解的 V、M₁、M₂ 即可做阶段二）
V(x, y, z) = 0.0
M₁(x, y, z) = 0.0
M₂(x, y, z) = 0.0

const to = TimerOutput()

gmsh.initialize()
@timeit to "open msh file" gmsh.open("/home/jason/vscode/mindlin_plate/msh/circular_tri3_3.msh")
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

    (∫wwdΩ => elements)(kʷʷ)
    (∫φwdΩ => elements)(kᵠʷ)
    (∫φφdΩ => elements)(kᵠᵠ)
    (∫κκdΩ => elements)(kᵠᵠ)

    (∫wqdΩ => elements)(fʷ)
    (∫φmdΩ => elements)(fᵠ)
end

@timeit to "assemble boundary penalties" begin
    bnd_names = sort([k for k in keys(entities) if startswith(k, "Γ")])
    isempty(bnd_names) && error("No boundary physical groups found (keys starting with 'Γ'). Available keys=$(collect(keys(entities)))")

    for name in bnd_names
        elements_Γ = getElements(nodes, entities[name])
        set𝝭!(elements_Γ)

        prescribe!(elements_Γ, :α => α, :g => w)
        (∫αwwdΓ => elements_Γ)(kʷʷ, fʷ)

        prescribe!(elements_Γ, :α => α, :g₁ => φ₁, :g₂ => φ₂, :n₁₁ => 1.0, :n₁₂ => 0.0, :n₂₂ => 1.0)
        (∫αφφdΓ => elements_Γ)(kᵠᵠ, fᵠ)

        # 自然边界项（阶段一为 0，但显式装配以保持结构一致）
        prescribe!(elements_Γ, :V => V)
        (∫wVdΓ => elements_Γ)(fʷ)

        prescribe!(elements_Γ, :M₁ => M₁, :M₂ => M₂)
        (∫φMdΓ => elements_Γ)(fᵠ)
    end
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
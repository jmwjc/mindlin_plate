
n = 10
npoly(n) = (n+1)*(n+2)÷2
nˢ = zeros(n,n)
Δʷ = zeros(n,n)
Δᵠ = zeros(n,n)
tf = Array{Bool}(undef,n,n)
for i in 1:n
    nʷ = npoly(i)
    for j in 1:n
        nᵠ = npoly(j)
        nᶠ = npoly(min(i,j+1))
        nˢ[i,j] = 2*nᵠ + nʷ - nᶠ 
        tf[i,j] = nˢ[i,j] >= nʷ
        Δʷ[i,j] = nˢ[i,j] - nʷ
        Δᵠ[i,j] = nˢ[i,j] - nᵠ
        max_nˢ = nˢ[i,j]
        min_nˢ = nʷ
        println("nʷ = $nʷ, nᵠ = $nᵠ, max_nˢ = $max_nˢ, min_nˢ = $min_nˢ")
    end
end

function print_max_min_nˢ(ndofʷ,ndofᵠ)
    i = ((1+8*ndofʷ)^0.5 - 3) ÷ 2
    nʷ = npoly(i)
    j = ((1+8*ndofᵠ)^0.5 - 3) ÷ 2
    nᵠ = npoly(j)
    nᶠ = npoly(min(i,j+1))
    max_nˢ = 2*ndofᵠ + ndofʷ - nᶠ 
    min_nˢ = nʷ
    println("nʷ = $nʷ, nᵠ = $nᵠ, nᶠ = $nᶠ, max_nˢ = $max_nˢ, min_nˢ = $min_nˢ")
end
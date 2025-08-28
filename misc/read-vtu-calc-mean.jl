using VTKDataIO, Statistics, LinearAlgebra, Plots

# --- 1. Ler o arquivo VTU ---
mesh = read_vtk("E:\\cobem2025\\results\\uasbcp1201.vtu")

points = mesh["points"]          # coordenadas [Npoints × 3]
cells = mesh["cells"]            # conectividade (lista de vetores de índices)

# --- 2. Variáveis (em cell_data) ---
U   = mesh["cell_data"]["uh"]        # velocidade (Ncélulas × 3, usar [:,1:2])
p   = mesh["cell_data"]["ph"]        # pressão
c   = mesh["cell_data"]["ch"]        # concentração
k   = mesh["cell_data"]["kh"]        # k turbulento
eps = mesh["cell_data"]["epsilon"]   # epsilon

# --- 3. Função para área de quadrilátero ---
function area_quad(coords)
    a, b, c_, d = coords
    A1 = 0.5 * norm(cross(b - a, c_ - a))
    A2 = 0.5 * norm(cross(d - a, c_ - a))
    return A1 + A2
end

# --- 4. Coordenada média e área de cada célula ---
yc = Float64[]
areas = Float64[]
for cell in cells
    coords = points[cell, :]
    push!(yc, mean(coords[:,2]))           # coordenada média em y
    push!(areas, area_quad(coords))        # área da célula
end

# --- 5. Discretizar eixo y em faixas ---
nbins = 50
ymin, ymax = extrema(yc)
ybins = range(ymin, ymax; length=nbins+1)

# --- 6. Calcular médias por faixa ---
Cmu = 0.09
ycenters = [(ybins[j] + ybins[j+1])/2 for j in 1:nbins]

Ux = fill(NaN, nbins)
Uy = fill(NaN, nbins)
pp = fill(NaN, nbins)
cc = fill(NaN, nbins)
kk = fill(NaN, nbins)
ee = fill(NaN, nbins)
nut = fill(NaN, nbins)

for j in 1:nbins
    ylow, yhigh = ybins[j], ybins[j+1]
    idx = findall((yc .>= ylow) .& (yc .< yhigh))
    if !isempty(idx)
        Aw = sum(areas[idx])
        Ux[j] = sum(U[idx,1] .* areas[idx]) / Aw
        Uy[j] = sum(U[idx,2] .* areas[idx]) / Aw
        pp[j] = sum(p[idx]     .* areas[idx]) / Aw
        cc[j] = sum(c[idx]     .* areas[idx]) / Aw
        kk[j] = sum(k[idx]     .* areas[idx]) / Aw
        ee[j] = sum(eps[idx]   .* areas[idx]) / Aw
        nut[j] = Cmu * kk[j]^2 / ee[j]
    end
end

# --- 7. Plotar perfis ---
plot(ycenters, Ux, label="Ux", xlabel="y", ylabel="Valor médio", lw=2)
plot!(ycenters, Uy, label="Uy", lw=2)
plot!(ycenters, pp, label="p", lw=2)
plot!(ycenters, cc, label="c", lw=2)
plot!(ycenters, kk, label="k", lw=2)
plot!(ycenters, ee, label="epsilon", lw=2)
plot!(ycenters, nut, label="νt", lw=2, linestyle=:dash)

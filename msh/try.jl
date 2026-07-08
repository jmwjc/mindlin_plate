using LinearAlgebra

"""
    generate_w_mesh(ndiv_q; output_dir=".")

给定 ndiv_q（对应 Q 场和 φ 场的 msh 文件 patchtest_tri3_<ndiv_q>.msh），
生成一个 w 场的 msh 文件，使得第一个秩条件 nʷ = n(n+1)/2 刚好满足（n_diff=0）。

输出文件: patchtest_tri3_w_<ndiv_q>.msh
"""
function generate_w_mesh(ndiv_q::Int; output_dir=".")
    # ── 1. 计算目标节点数 ──────────────────────────────────
    nˢ = (ndiv_q + 1)^2          # Q 场节点数
    n = floor(Int, 0.5 * ((1 + 8*nˢ)^0.5 - 3))  # 多项式阶数
    T_n = div(n * (n + 1), 2)    # 目标 w 场节点数 = n(n+1)/2

    println("ndiv_q = $ndiv_q, nˢ = $nˢ, n = $n, T_n = $T_n")

    # ── 2. 将 T_n 分解为 (nx+1)×(ny+1) = T_n ──────────────
    # 找最接近 sqrt(T_n) 的因子对作为 (nx+1) 和 (ny+1)
    best_nx1 = 1
    best_ny1 = T_n
    best_diff = T_n
    for nx1 in 1:T_n
        if T_n % nx1 == 0
            ny1 = T_n ÷ nx1
            diff = abs(nx1 - ny1)
            if diff < best_diff
                best_nx1 = nx1
                best_ny1 = ny1
                best_diff = diff
            end
        end
    end
    # nx1, ny1 是各方向的节点数
    nx1, ny1 = best_nx1, best_ny1
    if nx1 < ny1
        nx1, ny1 = ny1, nx1
    end
    # 划分数 = 节点数 - 1
    nx = nx1 - 1
    ny = ny1 - 1

    println("网格分解: $(nx1)×$(ny1) 个节点 (nx=$nx, ny=$ny)")
    @assert (nx+1)*(ny+1) == T_n "分解错误: ($(nx+1))*($(ny+1)) = $((nx+1)*(ny+1)) != $T_n"

    # ── 3. 生成节点坐标 ────────────────────────────────────
    # 单位正方形 [0,1]×[0,1]
    nodes = []
    node_id = 0
    # 角点
    push!(nodes, (0.0, 0.0))  # 1: 左下
    push!(nodes, (1.0, 0.0))  # 2: 右下
    push!(nodes, (1.0, 1.0))  # 3: 右上
    push!(nodes, (0.0, 1.0))  # 4: 左上
    node_id = 4

    # 底边 (不含角点)
    bottom_nodes = Int[]
    for i in 1:nx-1
        node_id += 1
        push!(nodes, (i/nx, 0.0))
        push!(bottom_nodes, node_id)
    end

    # 右边 (不含角点)
    right_nodes = Int[]
    for j in 1:ny-1
        node_id += 1
        push!(nodes, (1.0, j/ny))
        push!(right_nodes, node_id)
    end

    # 顶边 (不含角点, 从右到左)
    top_nodes = Int[]
    for i in nx-1:-1:1
        node_id += 1
        push!(nodes, (i/nx, 1.0))
        push!(top_nodes, node_id)
    end

    # 左边 (不含角点, 从上到下)
    left_nodes = Int[]
    for j in ny-1:-1:1
        node_id += 1
        push!(nodes, (0.0, j/ny))
        push!(left_nodes, node_id)
    end

    # 内部节点 (按行扫描)
    internal_nodes = []
    for j in 1:ny-1
        for i in 1:nx-1
            node_id += 1
            push!(nodes, (i/nx, j/ny))
            push!(internal_nodes, node_id)
        end
    end

    n_nodes = length(nodes)
    @assert n_nodes == T_n "节点数不匹配: $n_nodes != $T_n"

    # ── 4. 生成单元拓扑 ────────────────────────────────────
    # 每个矩形剖分为 2 个 Tri3
    # 节点编号映射: 结构化网格 (i,j) → 全局编号
    # 角点: 1(0,0), 2(nx,0), 3(nx,ny), 4(0,ny)
    # 底边: 5 ~ 4+(nx-1)
    # 右边: 5+(nx-1) ~ 4+(nx-1)+(ny-1)
    # 顶边: 5+(nx-1)+(ny-1) ~ 4+(nx-1)+(ny-1)+(nx-1)
    # 左边: 5+(nx-1)+(ny-1)+(nx-1) ~ 4+(nx-1)+(ny-1)+(nx-1)+(ny-1)
    # 内部: 5+(nx-1)+(ny-1)+(nx-1)+(ny-1) ~ end

    # 构建 (i,j) → id 的映射, i=0..nx, j=0..ny
    id_map = zeros(Int, nx+1, ny+1)
    id_map[1,1] = 1      # (0,0)
    id_map[nx+1,1] = 2   # (nx,0)
    id_map[nx+1,ny+1] = 3 # (nx,ny)
    id_map[1,ny+1] = 4   # (0,ny)

    # 底边 i=1..nx-1, j=0
    for (idx, i) in enumerate(1:nx-1)
        id_map[i+1, 1] = 4 + idx
    end
    # 右边 i=nx, j=1..ny-1
    offset = nx - 1
    for (idx, j) in enumerate(1:ny-1)
        id_map[nx+1, j+1] = 4 + offset + idx
    end
    # 顶边 i=nx-1..1, j=ny
    offset += ny - 1
    for (idx, i) in enumerate(nx-1:-1:1)
        id_map[i+1, ny+1] = 4 + offset + idx
    end
    # 左边 i=0, j=ny-1..1
    offset += nx - 1
    for (idx, j) in enumerate(ny-1:-1:1)
        id_map[1, j+1] = 4 + offset + idx
    end
    # 内部
    offset += ny - 1
    for j in 1:ny-1
        for i in 1:nx-1
            id_map[i+1, j+1] = 4 + offset + (j-1)*(nx-1) + i
        end
    end

    # 生成 Tri3 单元
    elements_tri = []
    for j in 0:ny-1
        for i in 0:nx-1
            n1 = id_map[i+1, j+1]
            n2 = id_map[i+2, j+1]
            n3 = id_map[i+2, j+2]
            n4 = id_map[i+1, j+2]
            # 两个三角形: (n1,n2,n3) 和 (n1,n3,n4)
            push!(elements_tri, (n1, n2, n3))
            push!(elements_tri, (n1, n3, n4))
        end
    end

    n_elements_tri = length(elements_tri)

    # ── 5. 边界单元 ────────────────────────────────────────
    # 底边 Γ¹: 节点 1 → 底边节点 → 2
    bottom_all = [1; bottom_nodes; 2]
    elements_Γ1 = [(bottom_all[i], bottom_all[i+1]) for i in 1:length(bottom_all)-1]

    # 右边 Γ²: 节点 2 → 右边节点 → 3
    right_all = [2; right_nodes; 3]
    elements_Γ2 = [(right_all[i], right_all[i+1]) for i in 1:length(right_all)-1]

    # 顶边 Γ³: 节点 3 → 顶边节点 → 4
    top_all = [3; top_nodes; 4]
    elements_Γ3 = [(top_all[i], top_all[i+1]) for i in 1:length(top_all)-1]

    # 左边 Γ⁴: 节点 4 → 左边节点 → 1
    left_all = [4; left_nodes; 1]
    elements_Γ4 = [(left_all[i], left_all[i+1]) for i in 1:length(left_all)-1]

    # Γ = Γ¹ ∪ Γ² ∪ Γ³ ∪ Γ⁴
    elements_Γ = vcat(elements_Γ1, elements_Γ2, elements_Γ3, elements_Γ4)

    n_Γ1 = length(elements_Γ1)
    n_Γ2 = length(elements_Γ2)
    n_Γ3 = length(elements_Γ3)
    n_Γ4 = length(elements_Γ4)
    n_Γ  = length(elements_Γ)
    n_Ω  = n_elements_tri

    # ── 6. 写入 Gmsh 4.1 格式 ──────────────────────────────
    filename = joinpath(output_dir, "patchtest_tri3_w_$(ndiv_q).msh")
    open(filename, "w") do io
        # 头部
        println(io, "\$MeshFormat")
        println(io, "4.1 0 8")
        println(io, "\$EndMeshFormat")

        # PhysicalNames
        println(io, "\$PhysicalNames")
        println(io, "6")
        println(io, "1 1 \"Γ¹\"")
        println(io, "1 2 \"Γ²\"")
        println(io, "1 3 \"Γ³\"")
        println(io, "1 4 \"Γ⁴\"")
        println(io, "1 6 \"Γ\"")
        println(io, "2 5 \"Ω\"")
        println(io, "\$EndPhysicalNames")

        # Entities
        println(io, "\$Entities")
        println(io, "4 5 1 0")
        # 点
        println(io, "1 0 0 0 0 ")
        println(io, "2 1 0 0 0 ")
        println(io, "3 1 1 0 0 ")
        println(io, "4 0 1 0 0 ")
        # 曲线
        println(io, "1 0 0 0 1 0 0 1 1 2 1 -2 ")
        println(io, "2 1 0 0 1 1 0 1 2 2 2 -3 ")
        println(io, "3 0 1 0 1 1 0 1 3 2 3 -4 ")
        println(io, "4 0 0 0 0 1 0 1 4 2 4 -1 ")
        println(io, "5 0 0 0 1 1 0 1 6 0 ")
        # 面
        println(io, "1 0 0 0 1 1 0 1 5 4 1 2 3 4 ")
        println(io, "\$EndEntities")

        # Nodes
        println(io, "\$Nodes")
        println(io, "10 $n_nodes 1 $n_nodes")
        # 角点
        println(io, "0 1 0 1")
        println(io, "1")
        println(io, "0 0 0")
        println(io, "0 2 0 1")
        println(io, "2")
        println(io, "1 0 0")
        println(io, "0 3 0 1")
        println(io, "3")
        println(io, "1 1 0")
        println(io, "0 4 0 1")
        println(io, "4")
        println(io, "0 1 0")

        # 底边
        if length(bottom_nodes) > 0
            println(io, "1 1 0 $(length(bottom_nodes))")
            for id in bottom_nodes
                println(io, id)
            end
            for id in bottom_nodes
                x, y = nodes[id]
                println(io, "$x $y 0")
            end
        end

        # 右边
        if length(right_nodes) > 0
            println(io, "1 2 0 $(length(right_nodes))")
            for id in right_nodes
                println(io, id)
            end
            for id in right_nodes
                x, y = nodes[id]
                println(io, "$x $y 0")
            end
        end

        # 顶边
        if length(top_nodes) > 0
            println(io, "1 3 0 $(length(top_nodes))")
            for id in top_nodes
                println(io, id)
            end
            for id in top_nodes
                x, y = nodes[id]
                println(io, "$x $y 0")
            end
        end

        # 左边
        if length(left_nodes) > 0
            println(io, "1 4 0 $(length(left_nodes))")
            for id in left_nodes
                println(io, id)
            end
            for id in left_nodes
                x, y = nodes[id]
                println(io, "$x $y 0")
            end
        end

        # 内部节点
        if length(internal_nodes) > 0
            println(io, "1 5 0 0")
            println(io, "2 1 0 $(length(internal_nodes))")
            for id in internal_nodes
                println(io, id)
            end
            for id in internal_nodes
                x, y = nodes[id]
                println(io, "$x $y 0")
            end
        end

        println(io, "\$EndNodes")

        # Elements
        total_elements = n_Γ1 + n_Γ2 + n_Γ3 + n_Γ4 + n_Γ + n_Ω
        println(io, "\$Elements")
        println(io, "6 $total_elements 1 $(total_elements)")

        # Γ¹ 边界线单元 (Seg2, element type 1)
        println(io, "1 1 1 $(n_Γ1)")
        for (idx, (n1, n2)) in enumerate(elements_Γ1)
            println(io, "$idx $n1 $n2")
        end

        # Γ² 边界线单元
        println(io, "1 2 1 $(n_Γ2)")
        for (idx, (n1, n2)) in enumerate(elements_Γ2)
            println(io, "$(n_Γ1+idx) $n1 $n2")
        end

        # Γ³ 边界线单元
        println(io, "1 3 1 $(n_Γ3)")
        for (idx, (n1, n2)) in enumerate(elements_Γ3)
            println(io, "$(n_Γ1+n_Γ2+idx) $n1 $n2")
        end

        # Γ⁴ 边界线单元
        println(io, "1 4 1 $(n_Γ4)")
        for (idx, (n1, n2)) in enumerate(elements_Γ4)
            println(io, "$(n_Γ1+n_Γ2+n_Γ3+idx) $n1 $n2")
        end

        # Γ 边界线单元 (所有边界) — 实体 tag=5, physicalTag=6
        println(io, "1 5 1 $(n_Γ)")
        for (idx, (n1, n2)) in enumerate(elements_Γ)
            println(io, "$(n_Γ1+n_Γ2+n_Γ3+n_Γ4+idx) 6 $n1 $n2")
        end

        # Ω 三角形单元 (Tri3, element type 2) — 实体 tag=1 (面), physicalTag=5
        println(io, "2 1 2 $(n_Ω)")
        for (idx, (n1, n2, n3)) in enumerate(elements_tri)
            println(io, "$(n_Γ1+n_Γ2+n_Γ3+n_Γ4+n_Γ+idx) $n1 $n2 $n3")
        end

        println(io, "\$EndElements")
    end

    println("已生成: $filename")
    println("  n_nodes = $n_nodes, n_elements = $n_elements_tri")
    println("  n_Γ1 = $n_Γ1, n_Γ2 = $n_Γ2, n_Γ3 = $n_Γ3, n_Γ4 = $n_Γ4")

    return filename
end

# ── 命令行入口 ──────────────────────────────────────────────
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        println("用法: julia try.jl <ndiv_q> [output_dir]")
        println("示例: julia try.jl 8")
        println("      julia try.jl 8 /path/to/output")
        exit(1)
    end
    ndiv_q = parse(Int, ARGS[1])
    output_dir = length(ARGS) >= 2 ? ARGS[2] : "."
    generate_w_mesh(ndiv_q, output_dir=output_dir)
end

import Gmsh: gmsh

function generateCircularMsh(filename::String; R::Float64=5.0, ndiv::Int=16)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2)

    a = R
    n_nodes = ndiv + 1

    # 创建点：圆心和轴上点
    gmsh.model.geo.addPoint(0.0, 0.0, 0.0, 1)  # 圆心
    gmsh.model.geo.addPoint(a, 0.0, 0.0, 2)     # x轴上点
    gmsh.model.geo.addPoint(0.0, a, 0.0, 3)     # y轴上点

    # 创建圆弧上的点：从x轴到y轴，均布角度
    arc_points = [2]  # 起始点
    for i in 1:ndiv
        angle = π/2 * i / ndiv
        x = a * cos(angle)
        y = a * sin(angle)
        point_id = 3 + i
        gmsh.model.geo.addPoint(x, y, 0.0, point_id)
        push!(arc_points, point_id)
    end
    push!(arc_points, 3)  # 结束点

    # 创建线
    gmsh.model.geo.addLine(1, 2, 1)  # Line 1: x轴 (0,0) -> (a,0)
    gmsh.model.geo.addLine(3, 1, 3)  # Line 3: y轴 (0,a) -> (0,0)

    # 创建圆弧线段
    arc_lines = []
    for i in 1:ndiv
        line_id = 3 + i
        gmsh.model.geo.addLine(arc_points[i], arc_points[i+1], line_id)
        push!(arc_lines, line_id)
    end

    # 创建曲线环和面
    curve_loop = vcat([1], arc_lines, [3])
    gmsh.model.geo.addCurveLoop(curve_loop, 1)
    gmsh.model.geo.addPlaneSurface([1], 1)

    gmsh.model.geo.synchronize()

    # 设置Transfinite
    gmsh.model.mesh.setTransfiniteCurve(1, n_nodes)
    for line in arc_lines
        gmsh.model.mesh.setTransfiniteCurve(line, n_nodes)
    end
    gmsh.model.mesh.setTransfiniteCurve(3, n_nodes)
    gmsh.model.mesh.setTransfiniteSurface(1)

    # 物理组
    # Ω: 域
    gmsh.model.addPhysicalGroup(2, [1], 1)
    gmsh.model.setPhysicalName(2, 1, "Ω")

    # Γ: 全部边界
    all_boundary = vcat([1], arc_lines, [3])
    gmsh.model.addPhysicalGroup(1, all_boundary, 2)
    gmsh.model.setPhysicalName(1, 2, "Γ")

    # Γ_circ: 圆弧边界
    gmsh.model.addPhysicalGroup(1, arc_lines, 3)
    gmsh.model.setPhysicalName(1, 3, "Γ_circ")

    # Γ_x: x轴边界
    gmsh.model.addPhysicalGroup(1, [1], 4)
    gmsh.model.setPhysicalName(1, 4, "Γ_x")

    # Γ_y: y轴边界
    gmsh.model.addPhysicalGroup(1, [3], 5)
    gmsh.model.setPhysicalName(1, 5, "Γ_y")

    # 生成网格
    gmsh.model.mesh.generate(2)

    # 保存
    gmsh.write(filename)
    gmsh.finalize()
    println("Generated: $filename (ndiv=$ndiv, arc segments=$ndiv)")
end

for n in 16:16
    generateCircularMsh("msh/circular_$n.msh", ndiv=n)
end

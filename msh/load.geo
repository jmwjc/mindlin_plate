Merge "patchtest_un_quad4_10.msh"; // 換成你實際的 msh 檔名
Delete Physical Groups;          // 擦除原本帶有內邊的錯誤物理組
CreateTopology;                  // 強制 Gmsh 擦除全場內邊，只提取最外圈邊界
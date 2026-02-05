a = 5.0;
n = 2;


Point(1) = {0.0, 0.0, 0.0};
Point(2) = {a/2, 0.0, 0.0};
Point(3) = {  a, 0.0, 0.0};
Point(4) = {a/2^0.5, a/2^0.5, 0.0};
Point(5) = {0.0,   a, 0.0};
Point(6) = {0.0, a/2, 0.0};
Point(7) = {a/2^0.5/2, a/2^0.5/2, 0.0};

Line(1) = {1,2};
Line(2) = {2,3};
Circle(3) = {3,1,4};
Circle(4) = {4,1,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line(7) = {7,2};
Line(8) = {7,6};
Line(9) = {7,4};

Curve Loop(1) = {1,-7,8,6};
Curve Loop(2) = {2,3,-9,7};
Curve Loop(3) = {9,4,5,-8};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Transfinite Curve{1,2,3,4,5,6,7,8,9} = n+1;
Transfinite Surface{1,2} ;
Transfinite Surface{3} ;

Physical Curve("Γᵇ") = {1,2};
Physical Curve("Γᵉ") = {3,4};
Physical Curve("Γˡ") = {5,6};
// 新增：整圈边界（合并所有边界曲线），用于程序中 entities["Γ"]
Physical Curve("Γ") = {1,2,3,4,5,6};
Physical Surface("Ω") = {1,2,3};
Physical Point("𝐴") = {1};

Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
//Mesh.SecondOrderIncomplete = 1;
Mesh 2;
//RecombineMesh;
//SetOrder 2;

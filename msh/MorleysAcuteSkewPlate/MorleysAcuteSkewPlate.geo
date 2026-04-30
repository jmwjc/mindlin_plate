
a = 100.0;
n = 7;
b = 3^0.5*a/2;
Point(1) = {  0.0, 0.0, 0.0};
Point(2) = {  a/2, 0.0, 0.0};
Point(3) = {    a, 0.0, 0.0};
Point(4) = {a+a/4, b/2, 0.0};
Point(5) = {a/2+a,   b, 0.0};
Point(6) = {    a,   b, 0.0};
Point(7) = {  a/2,   b, 0.0};
Point(8) = {  a/4, b/2, 0.0};
Point(9) = {3*a/4, b/2, 0.0};

Line(1)  = {1,2};
Line(2)  = {2,3};
Line(3)  = {3,4};
Line(4)  = {4,5};
Line(5)  = {5,6};
Line(6)  = {6,7};
Line(7)  = {7,8};
Line(8)  = {8,1};
Line(9)  = {2,9};
Line(10) = {9,4};
Line(11) = {9,6};
Line(12) = {9,8};

Curve Loop(1) = {1,9,12,8};
Curve Loop(2) = {2,3,-10,-9};
Curve Loop(3) = {10,4,5,-11};
Curve Loop(4) = {-12,11,6,7};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

Transfinite Curve{1,2,3,4,5,6,7,8,9,10,11,12} = n+1;
Transfinite Surface{1};
Transfinite Surface{2};
Transfinite Surface{3};
Transfinite Surface{4}Right;

Physical Curve("Γᵇ") = {1,2};
Physical Curve("Γʳ") = {3,4};
Physical Curve("Γᵗ") = {5,6};
Physical Curve("Γˡ") = {7,8};
Physical Surface("Ω") = {1,2,3,4};
Physical Point("𝐴") = {9};

Mesh.Algorithm = 1;
//Mesh.MshFileVersion = 2;
Mesh.SecondOrderIncomplete = 1;
Mesh 2;
RecombineMesh;
SetOrder 2;

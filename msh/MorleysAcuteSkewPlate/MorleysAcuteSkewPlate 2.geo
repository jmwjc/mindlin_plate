n = 32;
a = 100.0;
b = 3^0.5*a/2;
c = 100.0/n;
d = 3^0.5*a/2/n;

Point(1) = {  c*3/2,   d, 0.0};
Point(2) = {    a-c/2, d, 0.0};
Point(3) = {a/2+a-c*3/2,   b-d, 0.0};
Point(4) = {  a/2+c/2,   b-d, 0.0};

Line(1)  = {1,2};
Line(2)  = {2,3};
Line(3)  = {3,4};
Line(4)  = {4,1};

Curve Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

Transfinite Curve{1,2,3,4} = n-1;
Transfinite Surface{1};


Physical Curve("Γᵇ") = {1,2};
Physical Curve("Γʳ") = {3,4};
Physical Surface("Ω") = {1};
Mesh.Algorithm = 1;
//Mesh.MshFileVersion = 2;
Mesh.SecondOrderIncomplete = 1;
Mesh 2;
RecombineMesh;
SetOrder 2;

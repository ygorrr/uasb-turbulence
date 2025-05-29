lc = 0.1;
//+
Point(1) = {0, 0, 0, lc};
//+
Point(2) = {1, 0, 0, lc};
//+
Point(3) = {1, 1, 0, lc};
//+
Point(4) = {0, 1, 0, lc};
//+
Point(5) = {0.375, 0, 0, lc};
//+
Point(6) = {0.625, 0, 0, lc};
Point(7) = {0.5, 0, 0, lc};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {2, 6};
//+
Line(5) = {6, 7};
Line(6) = {7, 5};
//+
Line(7) = {5, 1};


Line Loop(1) = {7, 6, 5, 4, 3, 2, 1};
Plane Surface(1) = {1};

Physical Point("anchor") = {7};
Physical Line("inlet") = {5};
Physical Line("freeFlow") = {1, 2, 3, 4, 6};
//+
Physical Surface(7) = {1};

//Mesh.CharacteristicLengthMin = //+
MeshSize {5, 6, 7} = 0.01;

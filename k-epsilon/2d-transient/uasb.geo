h = 0.1;
L = 1.8;
H = 2.0;

Point(1) = {-L/2, 0, 0, h};
Point(2) = {L/2, 0, 0, h};
Point(3) = {L/2, H, 0, h};
Point(4) = {-L/2, H, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};

Surface(1) = {1};

Physical Line("bottom") = {1};
Physical Line("top") = {3};
Physical Line("sides") = {2, 4};

Mesh 2;

Save "uasb.msh";
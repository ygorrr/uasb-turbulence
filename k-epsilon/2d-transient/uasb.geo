h = 0.025;
L = 1.8;
H = 2*L;

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

Transfinite Line{1,3} = 1.5*50;
Transfinite Line{2,4} = (2)*50;
Transfinite Surface{1};
Recombine Surface{1};

Periodic Curve{4} = {2} Translate{L, 0, 0};

Physical Line("bottom") = {1};
Physical Line("top") = {3};
Physical Line("sides") = {2};
Physical Surface("domain") = {1};

Mesh 2;
Save "uasb.msh";
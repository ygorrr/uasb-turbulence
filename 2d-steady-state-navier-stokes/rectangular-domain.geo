// target mesh size (lc) close to the point
lc = 2e-1;

Point(1) = {0, 0, 0, lc};
Point(2) = {10, 0,  0, lc};
Point(3) = {10, 1, 0, lc};
Point(4) = {0,  1, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

// At this level, Gmsh knows everything to display the rectangular surface 1 and
// to mesh it. An optional step is needed if we want to group elementary
// geometrical entities into more meaningful groups, e.g. to define some
// mathematical ("domain", "boundary"), functional ("left wing", "fuselage") or
// material ("steel", "carbon") properties.
//
// Such groups are called "Physical Groups" in Gmsh. By default, if physical
// groups are defined, Gmsh will export in output files only mesh elements that
// belong to at least one physical group. (To force Gmsh to save all elements,
// whether they belong to physical groups or not, set `Mesh.SaveAll=1;', or
// specify `-save_all' on the command line.) Physical groups are also identified
// by tags, i.e. strictly positive integers, that should be unique per dimension
// (0D, 1D, 2D or 3D). Physical groups can also be given names.
//
// Here we define a physical curve that groups the left, bottom and right curves
// in a single group (with prescribed tag 5); and a physical surface with name
// "My surface" (with an automatic tag) containing the geometrical surface 1:

// Physical Curve(5) = {1, 2, 4};
// Physical Surface("My surface") = {1};

Physical Curve("bottom") = {1};
Physical Curve("right") = {2};
Physical Curve("top") = {3};
Physical Curve("left") = {4};

Physical Surface("domain") = {1};

// Now that the geometry is complete, you can
// - either open this file with Gmsh and select `2D' in the `Mesh' module to
//   create a mesh; then select `Save' to save it to disk in the default format
//   (or use `File->Export' to export in other formats);
// - or run `gmsh t1.geo -2` to mesh in batch mode on the command line.

// You could also uncomment the following lines in this script:

Mesh 2;
Save "rectangular-domain.msh";

// which would lead Gmsh to mesh and save the mesh every time the file is
// parsed. (To simply parse the file from the command line, you can use `gmsh
// t1.geo -')
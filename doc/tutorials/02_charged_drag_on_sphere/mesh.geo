SetFactory("OpenCASCADE");

upstream_distance = -0.045;
downstream_distance = 0.08;
domain_width = 0.035;

radius = 0.0095;

cs1 = 0.0003; // Cell size on the sphere
cs2 = 0.005;  // Cell size extending to domain boundaries

// ----- MESH ----- //
Point(1) = {upstream_distance, 0., 0., cs2};
Point(2) = {downstream_distance, 0., 0., cs2};
Point(3) = {downstream_distance, domain_width, 0., cs2};
Point(4) = {upstream_distance, domain_width, 0., cs2};

Point(5) = {-radius, 0., 0., cs1};
Point(6) = {0., radius, 0., cs1};
Point(7) = {radius, 0., 0., cs1};
Point(8) = {0., 0., 0., cs1};

Line(1) = {1,5};
Circle(2) = {5, 8, 6};
Circle(3) = {6, 8, 7};
Line(4) = {7,2};
Line(5) = {2,3};
Line(6) = {3,4};
Line(7) = {4,1};

Curve Loop(8) = {1,2,3,4,5,6,7};
Plane Surface(9) = {8};

Physical Curve("Sphere", 14) = {2, 3};
//+
Physical Curve("Inlet", 10) = {7};
//+
Physical Curve("Outlet", 11) = {5};
//+
Physical Curve("Top", 13) = {6};
//+
Physical Curve("Symmetry", 15) = {1,4};
//+
Physical Surface("Domain", 16) = {9};

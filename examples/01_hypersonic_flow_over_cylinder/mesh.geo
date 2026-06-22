cell_size = 0.005; // Set the cell size for the domain

Point(1) = {-1, 0, 0, cell_size};
Point(2) = {1, 0, 0, cell_size};
Point(3) = {1, 1, 0, cell_size};
Point(4) = {-1, 1, 0, cell_size};
Point(5) = {-0.1522, 0, 0, cell_size};
Point(6) = {0, 0.1522, 0, cell_size};
Point(7) = {0.1522, 0., 0, cell_size};
Point(8) = {0, 0, 0, cell_size};

Line(1) = {7, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {1, 5};

Circle(6) = {5, 8, 6};
Circle(7) = {6, 8, 7};

Curve Loop(1) = {3, 4, 5, 6, 7, 1, 2};
Plane Surface(1) = {1};

Physical Curve("Inlet", 8) = {4};
Physical Curve("Outlet", 9) = {2};
Physical Curve("Side", 10) = {3};
Physical Curve("Axis", 11) = {5, 1};
Physical Curve("Cylinder", 12) = {6, 7};
Physical Surface("Domain", 13) = {1};
cell_size = 1e-6;
domain_length = 1e-4;

Point(1) = {0, 0, 0, cell_size};
Point(2) = {domain_length, 0, 0, cell_size};
Line(1) = {1, 2};

Physical Curve("Domain", 2) = {1};
Physical Point("cathode", 3) = {1};
Physical Point("anode", 4) = {2};

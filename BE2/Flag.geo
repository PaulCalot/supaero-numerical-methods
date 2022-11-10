/*--------------------------*/
//    File Flag.geo
/*--------------------------*/

// Parameters
If (!Exists(h)) h=0.08; EndIf // Mesh size
If (!Exists(L)) L=2.; EndIf // Width
If (!Exists(l)) l=1.; EndIf // Height
If (!Exists(xc)) xc=L/3.; EndIf // x-center
If (!Exists(yc)) yc=l/2.; EndIf // y-center
If (!Exists(R)) R=l/6.; EndIf // Radius

// Points
Point(1) = {xc,yc,0,h};
Point(2) = {xc+R,yc,0,h};
Point(3) = {xc,yc+R,0,h};
Point(4) = {xc-R,yc,0,h};
Point(5) = {xc,yc-R,0,h};
Point(6) = {0,0,0,h};
Point(7) = {L,0,0,h};
Point(8) = {L,l,0,h};
Point(9) = {0,l,0,h};

// Curves
Circle(11) = {2,1,3};
Circle(12) = {3,1,4};
Circle(13) = {4,1,5};
Circle(14) = {5,1,2};
Line(15) = {6,7};
Line(16) = {7,8};
Line(17) = {8,9};
Line(18) = {9,6};

// Paths
Line Loop(21) = {11,12,13,14};
Line Loop(22) = {15,16,17,18};

// Surfaces
Plane Surface(31) = {-21,22};

// Physical Groups
Physical Line("Bottom") = {15};
Physical Line("Right") = {16};
Physical Line("Top") = {17};
Physical Line("Left") = {18};
Physical Line("Hole") = {11,12,13,14};
Physical Surface("Omega") = {31};


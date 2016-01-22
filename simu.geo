h=0.1;

Function hole
p0 = newp; Point(p0) = {x,y,0,h};
p1 = newp; Point(p1) = {x+r,y,0,h};
p2 = newp; Point(p2) = {x,y+r,0,h};
p3 = newp; Point(p3) = {x-r,y,0,h};
p4 = newp; Point(p4) = {x,y-r,0,h};

r1 = newreg; Circle(r1) = {p1,p0,p2};
r2 = newreg; Circle(r2) = {p2,p0,p3};
r3 = newreg; Circle(r3) = {p3,p0,p4};
r4 = newreg; Circle(r4) = {p4,p0,p1};

A = newreg; Line Loop(A) = {r1,r2,r3,r4};
Return


xmin=-1;
ymin=-1;
xmax=1;
ymax=3;
r=0.05;
Point(1) = {xmin,ymin, 0, h};
Point(2) = {xmax,ymin, 0, h};
Point(3) = {xmax,ymax, 0, h};
Point(4) = {xmin,ymax, 0, h};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(6) = {4, 1, 2, 3};

/// Generate holes
For Ax In {xmin+(xmax-xmin)/6.:xmin+5*(xmax-xmin)/6.:6*r}
For Ay In {ymin+(ymax-ymin)/6.:ymin+5*(ymax-ymin)/6.:6*r}
x = Ax+(Rand(2*r)-r);
y = Ay+(Rand(2*r)-r);
Call hole;
tabL[]+=r1;
tabL[]+=r2;
tabL[]+=r3;
tabL[]+=r4;
tab[]+=A;
EndFor
EndFor

For i In {0:10}
Printf("%g - %g", i, tab[i]);
EndFor


Plane Surface(6) = {6,tab[]};

Physical Surface("omega") = {6};
Physical Point("p") = {1};
Physical Line("g") = {3};
Physical Line("null") = {2,4,tabL[]};

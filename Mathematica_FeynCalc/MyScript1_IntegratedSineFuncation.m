Lx=10;
Kmin=-40;
Kmax=40;
PreFactor=0.02;


Plot3D[ (1 - Cos[PreFactor*(x*x + Lx*Lx - 2.0*x*Lx + y*y)]) / ((x*x + Lx*Lx - 2.0*x*Lx + y*y)^2 ) ,  {x,Kmin,Kmax}, {y,Kmin,Kmax} ]

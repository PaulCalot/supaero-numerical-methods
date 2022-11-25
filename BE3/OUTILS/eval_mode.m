function [E,H] = eval_mode(t,x,y,m,n,c0)

mu0 = 4*pi*1E-7;
eps0 = 1/mu0/c0^2;

w = c0*sqrt((pi*m)^2+(pi*n)^2);

E = pi/w/eps0*sin(w*t)*[n*cos(m*pi*x)*sin(n*pi*y);-m*sin(m*pi*x)*cos(n*pi*y) ];
H = -cos(w*t)*cos(m*pi*x)*cos(n*pi*y);


end 
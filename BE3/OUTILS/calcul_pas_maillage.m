function [hmin, hmax] = calcul_pas_maillage(mesh)

hmin = 1E9;
hmax = 0;

for T = 1:mesh.Nbtri
 xT = mesh.coor(1,mesh.Tri(1:3,T));
 yT = mesh.coor(2,mesh.Tri(1:3,T));
 AireT = Aire(xT,yT);
 hmin = min(hmin,AireT);
 hmax = max(hmax,AireT);
end  
hmin = sqrt(hmin);
hmax = sqrt(hmax);
end 
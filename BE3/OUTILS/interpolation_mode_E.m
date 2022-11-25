function  E = interpolation_mode_E(t,m,n,c0,mesh,ordre,N,ri,si)

% (ri,si) les N points d'interpolation sur le triangle de référence


for T = 1:mesh.Nbtri
 xT = mesh.coor(1,mesh.Tri(1:3,T));
 yT = mesh.coor(2,mesh.Tri(1:3,T));
 ind = N*(T-1);
 for l = 1:N
  [x,y] = eval_FT(ri(l),si(l),xT,yT);
  [Ep,Hp] = eval_mode(t,x,y,m,n,c0);
  E(2*ind+2*(l-1)+1,1) = Ep(1);
  E(2*ind+2*(l-1)+2,1) = Ep(2);
  
 end 
 
end
end 
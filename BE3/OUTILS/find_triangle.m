function [Ts, xref, yref ] = find_triangle(mesh,xs,ys)

Ts =1; 
trouve = 0;
while trouve == 0
 xT = mesh.coor(1,mesh.Tri(1:3,Ts));
 yT = mesh.coor(2,mesh.Tri(1:3,Ts));
 [xref,yref] = eval_inv_FT(xs,ys,xT,yT);
 if xref>=0 && xref<=1 && yref >=0 && yref <=1 && xref+yref-1<=0
  trouve =1
 else
  Ts = Ts+1; 
 end 
end
end 
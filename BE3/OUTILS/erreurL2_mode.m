function [erreurL2_E,erreurL2_H] = erreurL2_mode(E,H,mx,my,c0,niter,dt,A,ordre,mesh)


N = (ordre +1)*(ordre+2)/2;

tE = niter*dt;
tH= tE + dt/2;

% Choix d'une formule d'intégration numérique
ordreq = 4*ordre;
Nq = twb_rule_n (ordreq);
if (Nq ==-1)  %% ordre de quadrature inexistant
 ordreq = ordreq + 1;
 Nq = twb_rule_n (ordreq);
end  
xq = twb_rule_x ( ordreq );
yq = twb_rule_y ( ordreq );
wq = twb_rule_w ( ordreq );

erreurL2_E = 0;
normL2_E_ex = 0;
erreurL2_H = 0;
normL2_H_ex = 0;
for T = 1:mesh.Nbtri
 xT = mesh.coor(1,mesh.Tri(1:3,T));
 yT = mesh.coor(2,mesh.Tri(1:3,T));
 DF = DFT(xT,yT);
 JAC = abs(det(DF));
 ET = E(2*N*(T-1)+1:2*N*(T-1)+2*N,1);
 HT = H(N*(T-1)+1:N*(T-1)+N,1);
 for n = 1:Nq
   [x,y] = eval_FT(xq(n),yq(n),xT,yT); 
   [Eex,Hex] = eval_mode(tE,x,y,mx,my,c0);
   [Ep,Hp] = Eval_champs_ponctuels(ET,HT,xq(n),yq(n),ordre,A);
   erreurL2_E = erreurL2_E + JAC*wq(n)*((Ep-Eex)'*(Ep-Eex));
   normL2_E_ex = normL2_E_ex + JAC*wq(n)*(Eex'*Eex);
  
   [Eex,Hex] = eval_mode(tH,x,y,mx,my,c0);
   erreurL2_H = erreurL2_H + JAC*wq(n)*((Hp-Hex)^2);
   normL2_H_ex = normL2_H_ex + JAC*wq(n)*(Hex^2);
 end 
 
 end 
% normL2_E_ex;
% normL2_H_ex;
% erreurL2_E;
% erreurL2_H;
 erreurL2_E = sqrt(erreurL2_E/normL2_E_ex);
 erreurL2_H = sqrt(erreurL2_H/normL2_H_ex);
 


end 

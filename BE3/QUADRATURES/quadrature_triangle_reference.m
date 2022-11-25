function [Nq,xq,yq,wq] = quadrature_triangle_reference(ordreq)


Nq = twb_rule_n (ordreq);
if (Nq ==-1)  %% ordre de quadrature inexistant
 ordreq = ordreq + 1;
 Nq = twb_rule_n (ordreq);
end  
xq = twb_rule_x ( ordreq );
yq = twb_rule_y ( ordreq );
wq = twb_rule_w ( ordreq );
end 

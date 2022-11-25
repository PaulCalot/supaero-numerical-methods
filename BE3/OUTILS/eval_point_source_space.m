function Js = eval_point_source_space(x,y,xc,yc,r0)

%% (xc,yc) = centre du point source
%% sigma = ecart type

r = sqrt((x-xc)^2+(y-yc)^2);
Js = exp(-r^2/r0^2);

end 
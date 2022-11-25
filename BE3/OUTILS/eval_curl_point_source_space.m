function curl_Js = eval_curl_point_source_space(x,y,xc,yc,r0)

%% (xc,yc) = centre du point source
%% r0 = "rayon" de la source

r = sqrt((x-xc)^2+(y-yc)^2);
Js = exp(-r^2/r0^2);
curl_Js = 2*Js/r0^2*[-(y-yc); (x-xc)];

end 
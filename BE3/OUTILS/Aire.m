function  AireT = Aire(xT,yT)

V1 = [xT(2)-xT(1);yT(2)-yT(1)];
V2 = [xT(3)-xT(1);yT(3)-yT(1)];

AireT = abs(V1(1)*V2(2)-V1(2)*V2(1))/2;

end 
p1=0.5;
p2=0.0;
p3=1-p1-p2;
gG=1;
gB=1;

pi=0:0.0001:1;

Ab=p3.*(1-pi);
Ag=p1+(p2+p3).*pi;

AbP=-p3;
AgP=p2+p3;

X=gG*gB*(AbP*Ab.^2+AgP*Ag.^2+Ag.*Ab.*(AgP+AbP))...+(gG^2+gB^2)*(AbP*Ag.^2+AgP*Ab.^2);

plot(pi,X)

pi=0:0.0001:1;
lp=length(pi);


p1=0.5;
p2=0.3;
p3=0.2;
gG=0.9;
gB=0.3;
y=2;
D=1;
r=1.01;

[w,l]=meshgrid(0.5:0.001:1,0:0.001:1);
phatAH=-(w.*l.*p1*(y*gG-D*r)+(1-w).*(1-l).*p3*(y*gB-D*r))./(w.*l.*(p2+p3)*(y*gG-D*r)+(1-w).*(1-l).*(-p3)*(y*gB-D*r));
phatCH=-(w.*l.*p3*(y*gG-D*r)+(1-w).*(1-l).*p1*(y*gB-D*r))./(w.*l.*(-p3)*(y*gG-D*r)+(1-w).*(1-l).*(p2+p3)*(y*gB-D*r));
phatAL=-((1-w).*l.*p1*(y*gG-D*r)+w.*(1-l).*p3*(y*gB-D*r))./((1-w).*l.*(p2+p3)*(y*gG-D*r)+w.*(1-l).*(-p3)*(y*gB-D*r));
phatCL=-((1-w).*l.*p3*(y*gG-D*r)+w.*(1-l).*p1*(y*gB-D*r))./((1-w).*l.*(-p3)*(y*gG-D*r)+w.*(1-l).*(p2+p3)*(y*gB-D*r));

[c,h]=contour(w,l,phatCL,-0.8:0.2:0.8);
clabel(c,h);

w=0.7;
l=0.6;

wbarBH=-(1-l)*(y*gB-D*r)/(l*(y*gG-D*r)-(1-l)*(y*gB-D*r));
wbarBL=l*(y*gG-D*r)/(l*(y*gG-D*r)-(1-l)*(y*gB-D*r));
for i=1:lp
    wbarAH(i)=-(1-l)*p3*(1-pi(i))*(y*gB-D*r)/(l*(p1+(p2+p3)*pi(i))*(y*gG-D*r)-(1-l)*(p1+(p2+p3)*pi(i))*(y*gB-D*r));
    wbarCH(i)=-(1-l)*(p1+(p2+p3)*pi(i))*(y*gB-D*r)/(l*p3*(1-pi(i))*(y*gG-D*r)-(1-l)*(p1+(p2+p3)*pi(i))*(y*gB-D*r));
    wbarAL(i)=l*(p1+(p2+p3)*pi(i))*(y*gG-D*r)/(l*(p1+(p2+p3)*pi(i))*(y*gG-D*r)-(1-l)*p3*(1-pi(i))*(y*gB-D*r));
    wbarCL(i)=l*p3*(1-pi(i))*(y*gG-D*r)/(l*p3*(1-pi(i))*(y*gG-D*r)-(1-l)*(p1+(p2+p3)*pi(i))*(y*gB-D*r));
end

plot(pi,wbarAH,pi,wbarCH,pi,wbarAL,pi,wbarCL)
legend('AH','CH','AL','CL')
line([0 1],[wbarBH wbarBH],'LineStyle','-.')
line([0 1],[wbarBL wbarBL],'LineStyle','--')

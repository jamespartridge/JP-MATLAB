function [piHeq,piLeq,Rah,Rbh,Rch,Ral,Rbl,Rcl]=FP(w,par)

gG=par(1);         ...prob. pays off if good
gB=par(2);         ...prob. pays off if bad
y=par(3);          ...output
D=par(4);          ...investment size
r=par(5);          ...risk free rate
g1=par(6);         ...Pr(A|G)=g1+(g2+g3)pi
g2=par(7);         ...Pr(B|G)=g2(1-pi)
g3=par(8);         ...Pr(C|G)=g3(1-pi)
b1=par(9);         ...Pr(A|L)=b1+(b2+b3)pi
b2=par(10);        ...Pr(B|L)=b2(1-pi)
b3=par(11);        ...Pr(C|L)=b3(1-pi)
alf=par(12);       ...c(pi)=1/alpha * pi^alpha
l=par(13);

%pi grid
pi=0:0.00001:1;
lp=length(pi);

% interest rates
R=zeros(3,2,lp);
R(1,1,:)=r*(w*l*(g1+(g2+g3).*pi)+(1-w)*(1-l)*(b2+b3).*pi)./(w*l*(g1+(g2+g3).*pi)*gG+(1-w)*(1-l)*(b2+b3).*pi*gB);
R(2,1,:)=r*(w*l*g2.*(1-pi)+(1-w)*(1-l)*b2.*(1-pi))./(w*l*g2.*(1-pi)*gG+(1-w)*(1-l)*b2.*(1-pi)*gB);
R(3,1,:)=r*(w*l*g3.*(1-pi)+(1-w)*(1-l)*(b1+b3.*(1-pi)))./(w*l*g3.*(1-pi)*gG+(1-w)*(1-l)*(b1+b3.*(1-pi))*gB);
R(1,2,:)=r*((1-w)*l*(g1+(g2+g3).*pi)+w*(1-l)*(b2+b3).*pi)./((1-w)*l*(g1+(g2+g3).*pi)*gG+w*(1-l)*(b2+b3).*pi*gB);
R(2,2,:)=r*((1-w)*l*g2.*(1-pi)+w*(1-l)*b2.*(1-pi))./((1-w)*l*g2.*(1-pi)*gG+w*(1-l)*b2.*(1-pi)*gB);
R(3,2,:)=r*((1-w)*l*g3.*(1-pi)+w*(1-l)*(b1+b3.*(1-pi)))./((1-w)*l*g3.*(1-pi)*gG+w*(1-l)*(b1+b3.*(1-pi))*gB);



%individual terms of marginal benefit function
Hx(1,:)=(y-D*R(1,1,:)>0).*(y-D*R(1,1,:)).*((w*l*(g2+g3)*gG+(1-w)*(1-l)*(b2+b3)*gB)/(w*l+(1-w)*(1-l)));
Hx(2,:)=(y-D*R(2,1,:)>0).*(y-D*R(2,1,:)).*((w*l*(-g2)*gG+(1-w)*(1-l)*(-b2)*gB)/(w*l+(1-w)*(1-l)));
Hx(3,:)=(y-D*R(3,1,:)>0).*(y-D*R(3,1,:)).*((w*l*(-g3)*gG+(1-w)*(1-l)*(-b3)*gB)/(w*l+(1-w)*(1-l)));
Lx(1,:)=(y-D*R(1,2,:)>0).*(y-D*R(1,2,:)).*(((1-w)*l*(g2+g3)*gG+w*(1-l)*(b2+b3)*gB)/((1-w)*l+w*(1-l)));
Lx(2,:)=(y-D*R(2,2,:)>0).*(y-D*R(2,2,:)).*(((1-w)*l*(-g2)*gG+w*(1-l)*(-b2)*gB)/((1-w)*l+w*(1-l)));
Lx(3,:)=(y-D*R(3,2,:)>0).*(y-D*R(3,2,:)).*(((1-w)*l*(-g3)*gG+w*(1-l)*(-b3)*gB)/((1-w)*l+w*(1-l)));

%marginal benefit function
xH=Hx(1,:)+Hx(2,:)+Hx(3,:);
xL=Lx(1,:)+Lx(2,:)+Lx(3,:);
% firm's optimal pi
piH=xH.^(1/(alf-1));
piL=xL.^(1/(alf-1));
% enforce constraints
piH(xH<0)=0; piH(xH>1)=1;
piL(xL<0)=0; piL(xL>1)=1;
% eq is where investor belief = firm choice
[~,eqH]=min(abs(piH-pi));
[~,eqL]=min(abs(piL-pi));
% eq pi
piHeq=pi(eqH); piLeq=pi(eqL);
% interest rates at eq
Rah=R(1,1,eqH); Rbh=R(2,1,eqH); Rch=R(3,1,eqH);
Ral=R(1,2,eqL); Rbl=R(2,2,eqL); Rcl=R(3,2,eqL);

return    
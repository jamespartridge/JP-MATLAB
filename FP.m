function [piHeq,piLeq,Rah,Rbh,Rch,Ral,Rbl,Rcl]=FP(w,par)

gG=par(1);         ...prob. pays off if good
gB=par(2);         ...prob. pays off if bad
y=par(3);             ...output
D=par(4);             ...investment size
r=par(5);          ...risk free rate
g1=par(6);          ...Pr(A|G)=g1+(g2+g3)pi
g2=par(7);          ...Pr(B|G)=g2(1-pi)
g3=par(8);          ...Pr(C|G)=g3(1-pi)
b1=par(9);          ...Pr(A|L)=b1+(b2+b3)pi
b2=par(10);          ...Pr(B|L)=b2(1-pi)
b3=par(11);          ...Pr(C|L)=b3(1-pi)
alf=par(12);           ...c(pi)=1/alpha * pi^alpha
l=par(13);

pi=0:0.0001:1;
lp=length(pi);

%omega cutoffs and interest rates
wbar=zeros(3,2,lp);
wbar(2,1,:)=-(1-l)*(y*gB-D*r)/(l*(y*gG-D*r)-(1-l)*(y*gB-D*r));
wbar(2,2,:)=l*(y*gG-D*r)/(l*(y*gG-D*r)-(1-l)*(y*gB-D*r));
R=zeros(3,2,lp);
R(2,1,:)=r*(w*l+(1-w)*(1-l))/(w*l*gG+(1-w)*(1-l)*gB);
R(2,2,:)=r*((1-w)*l+w*(1-l))/((1-w)*l*gG+w*(1-l)*gB);
for i=1:lp
    wbar(1,1,i)=-(1-l)*b3*(1-pi(i))*(y*gB-D*r)/(l*(g1+(g2+g3)*pi(i))*(y*gG-D*r)-(1-l)*(b1+(b2+b3)*pi(i))*(y*gB-D*r));
    wbar(3,1,i)=-(1-l)*(b1+(b2+b3)*pi(i))*(y*gB-D*r)/(l*g3*(1-pi(i))*(y*gG-D*r)-(1-l)*(b1+(b2+b3)*pi(i))*(y*gB-D*r));
    wbar(1,2,i)=l*(g1+(g2+g3)*pi(i))*(y*gG-D*r)/(l*(g1+(g2+g3)*pi(i))*(y*gG-D*r)-(1-l)*b3*(1-pi(i))*(y*gB-D*r));
    wbar(3,2,i)=l*g3*(1-pi(i))*(y*gG-D*r)/(l*g3*(1-pi(i))*(y*gG-D*r)-(1-l)*(b1+(b2+b3)*pi(i))*(y*gB-D*r));
    R(1,1,i)=r*(w*l*(g1+(g2+g3)*pi(i))+(1-w)*(1-l)*b3*(1-pi(i)))/(w*l*(g1+(g2+g3)*pi(i))*gG+(1-w)*(1-l)*b3*(1-pi(i))*gB);
    R(3,1,i)=r*(w*l*g3*(1-pi(i))+(1-w)*(1-l)*(b1+(b2+b3)*pi(i)))/(w*l*g3*(1-pi(i))*gG+(1-w)*(1-l)*(b1+(b2+b3)*pi(i))*gB);
    R(1,2,i)=r*((1-w)*l*(g1+(g2+g3)*pi(i))+w*(1-l)*b3*(1-pi(i)))/((1-w)*l*(g1+(g2+g3)*pi(i))*gG+w*(1-l)*b3*(1-pi(i))*gB);
    R(3,2,i)=r*((1-w)*l*b3*(1-pi(i))+w*(1-l)*(g1+(g2+g3)*pi(i)))/((1-w)*l*b3*(1-pi(i))*gG+w*(1-l)*(g1+(g2+g3)*pi(i))*gB);
end

%pi cutoffs
phat=zeros(3,2);
phat(1,1)=-(w*l*g1*(y*gG-D*r)+(1-w)*(1-l)*b3*(y*gB-D*r))/(w*l*(g2+g3)*(y*gG-D*r)+(1-w)*(1-l)*(-b3)*(y*gB-D*r));
phat(3,1)=-(w*l*g3*(y*gG-D*r)+(1-w)*(1-l)*b1*(y*gB-D*r))/(w*l*(-g3)*(y*gG-D*r)+(1-w)*(1-l)*(b2+b3)*(y*gB-D*r));
phat(1,2)=-((1-w)*l*g1*(y*gG-D*r)+w*(1-l)*b3*(y*gB-D*r))/((1-w)*l*(g2+g3)*(y*gG-D*r)+w*(1-l)*(-b3)*(y*gB-D*r));
phat(3,2)=-((1-w)*l*g3*(y*gG-D*r)+w*(1-l)*b1*(y*gB-D*r))/((1-w)*l*(-g3)*(y*gG-D*r)+w*(1-l)*(b2+b3)*(y*gB-D*r));


%individual terms of marginal benefit function
Hx=zeros(3,lp);
Lx=zeros(3,lp);
for i=1:lp
    Hx(1,i)=((w*l*(g2+g3)*gG+(1-w)*(1-l)*(-b3)*gB)/(w*l+(1-w)*(1-l)))*(y-D*R(1,1,i));
    Hx(2,i)=((w*l*(-g2)*gG+(1-w)*(1-l)*(-b2)*gB)/(w*l+(1-w)*(1-l)))*(y-D*R(2,1,1));
    Hx(3,i)=((w*l*(-g3)*gG+(1-w)*(1-l)*(b2+b3)*gB)/(w*l+(1-w)*(1-l)))*(y-D*R(3,1,i));
    Lx(1,i)=(((1-w)*l*(g2+g3)*gG+w*(1-l)*(-b3)*gB)/((1-w)*l+w*(1-l)))*(y-D*R(1,2,i));
    Lx(2,i)=(((1-w)*l*(-g2)*gG+w*(1-l)*(-b2)*gB)/((1-w)*l+w*(1-l)))*(y-D*R(2,2,1));
    Lx(3,i)=(((1-w)*l*(-g3)*gG+w*(1-l)*(b2+b3)*gB)/((1-w)*l+w*(1-l)))*(y-D*R(3,2,i));
end

%check conditions to construct marginal benefit function
%omega > cutoff for H, omega < cutoff for L
%pi > cutoff for A, pi < cutoff for C
xH=zeros(1,lp);
xL=zeros(1,lp);
% for i=1:lp
%     %High
%     if w>wbar(1,1,i) && pi(i)>phat(1,1),    ...A?
%         xH(i)=Hx(1,i);  end
%     if w>wbar(2,1,i),                       ...B?
%         xH(i)=xH(i)+Hx(2,i); end
%     if w>wbar(3,1,i) && pi(i)<phat(3,1),    ...C?
%         xH(i)=xH(i)+Hx(3,i); end
%     %Low
%     if w<wbar(1,2,i) && pi(i)>phat(1,2),    ...A?
%         xL(i)=Lx(1,i);  end
%     if w<wbar(2,1,i),                       ...B?
%         xL(i)=xL(i)+Lx(2,i); end
%     if w<wbar(3,2,i) && pi(i)<phat(3,2),    ...C?
%         xL(i)=xL(i)+Lx(3,i); end
% end

for i=1:lp
    %High
    if y-D*R(1,1,i)>0,    ...A?
        xH(i)=Hx(1,i);  end
    if y-D*R(2,1,i)>0,    ...B?
        xH(i)=xH(i)+Hx(2,i); end
    if y-D*R(3,1,i)>0,    ...C?
        xH(i)=xH(i)+Hx(3,i); end
    %Low
    if y-D*R(1,2,i)>0,    ...A?
        xL(i)=Lx(1,i);  end
    if y-D*R(2,2,i)>0,    ...B?
        xL(i)=xL(i)+Lx(2,i); end
    if y-D*R(3,2,i)>0,    ...C?
        xL(i)=xL(i)+Lx(3,i); end
end


%firm's optimal pi
piH=zeros(1,lp);
piL=zeros(1,lp);
for i=1:lp
%     piH(i)=1-xH(i)^(-0.5);
%     piL(i)=1-xL(i)^(-0.5);
    piH(i)=xH(i)^(1/(alf-1));
    piL(i)=xL(i)^(1/(alf-1));
    if xH(i)<0, piH(i)=0; end
    if xL(i)<0, piL(i)=0; end
    if xH(i)>1, piH(i)=1; end
    if xL(i)>1, piL(i)=1; end
end
[~,eqH]=min(abs(piH-pi));
[~,eqL]=min(abs(piL-pi));

Rah=R(1,1,eqH); Rbh=R(2,1,eqH); Rch=R(3,1,eqH);
Ral=R(1,2,eqL); Rbl=R(2,2,eqL); Rcl=R(3,2,eqL);

piHeq=pi(eqH); piLeq=pi(eqL);

return

% MBHbyw=(y-D*R(1,1,eqH))*l*(1-l)*(gG*(g2+g3)-gB*(-b3))/((w*l+(1-w)*(1-l))^2)...
%         +((w*l*gG*(g2+g3)+(1-w)*(1-l)*gB*(-b3))/(w*l+(1-w)*(1-l)))*(-D*r*l*(1-l)*(g1+(g2+g3)*pi(eqH))*b3*(1-pi(eqH))*(gB-gG)/((w*l*(g1+(g2+g3)*pi(eqH))*gG+(1-w)*(1-l)*b3*(1-pi(eqH))*gB)^2))...
%       +(y-D*R(1,1,eqH))*l*(1-l)*(gG*(-g2)-gB*(-b2))/((w*l+(1-w)*(1-l))^2)...
%         +((w*l*gG*(-g2)+(1-w)*(1-l)*gB*(-b2))/(w*l+(1-w)*(1-l)))*(-D*r*l*(1-l)*(gB-gG)/((w*l*gG+(1-w)*(1-l)*gB)^2))...
%       +(y-D*R(1,1,eqH))*l*(1-l)*(gG*(-g3)-gB*(b2+b3))/((w*l+(1-w)*(1-l))^2)...
%         +((w*l*gG*(-g3)+(1-w)*(1-l)*gB*(b2+b3))/(w*l+(1-w)*(1-l)))*(-D*r*l*(1-l)*g3*(1-pi(eqH))*(b1+(b2+b3)*pi(eqH))*(gB-gG)/((w*l*b3*(1-pi(eqH))*gG+(1-w)*(1-l)*(b1+(b2+b3)*pi(eqH))*gB)^2));
    
    
    
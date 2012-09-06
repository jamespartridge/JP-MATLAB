pi=0:0.0001:1;
lp=length(pi);
x=zeros(1,lp);

w=0.7;
l=0.5;
p1=0.5;
p2=0.3;
p3=0.2;
gG=0.9;
gB=0.;
y=1.5;
D=1;
r=1.01;

%omega cutoffs and interest rates
wbar=zeros(3,2,lp);
wbar(2,1,:)=-(1-l)*(y*gB-D*r)/(l*(y*gG-D*r)-(1-l)*(y*gB-D*r));
wbar(2,2,:)=l*(y*gG-D*r)/(l*(y*gG-D*r)-(1-l)*(y*gB-D*r));
R=zeros(3,2,lp);
R(2,1,:)=r*(w*l+(1-w)*(1-l))/(w*l*gG+(1-w)*(1-l)*gB);
R(2,2,:)=R(2,1);
for i=1:lp
    wbar(1,1,i)=-(1-l)*p3*(1-pi(i))*(y*gB-D*r)/(l*(p1+(p2+p3)*pi(i))*(y*gG-D*r)-(1-l)*(p1+(p2+p3)*pi(i))*(y*gB-D*r));
    wbar(3,1,i)=-(1-l)*(p1+(p2+p3)*pi(i))*(y*gB-D*r)/(l*p3*(1-pi(i))*(y*gG-D*r)-(1-l)*(p1+(p2+p3)*pi(i))*(y*gB-D*r));
    wbar(1,2,i)=l*(p1+(p2+p3)*pi(i))*(y*gG-D*r)/(l*(p1+(p2+p3)*pi(i))*(y*gG-D*r)-(1-l)*p3*(1-pi(i))*(y*gB-D*r));
    wbar(3,2,i)=l*p3*(1-pi(i))*(y*gG-D*r)/(l*p3*(1-pi(i))*(y*gG-D*r)-(1-l)*(p1+(p2+p3)*pi(i))*(y*gB-D*r));
    R(1,1,i)=r*(w*l*(p1+(p2+p3)*pi(i))+(1-w)*(1-l)*p3*(1-pi(i)))/(w*l*(p1+(p2+p3)*pi(i))*gG+(1-w)*(1-l)*p3*(1-pi(i))*gB);
    R(3,1,i)=r*(w*l*p3*(1-pi(i))+(1-w)*(1-l)*(p1+(p2+p3)*pi(i)))/(w*l*p3*(1-pi(i))*gG+(1-w)*(1-l)*(p1+(p2+p3)*pi(i))*gB);
    R(1,2,i)=r*((1-w)*l*(p1+(p2+p3)*pi(i))+w*(1-l)*p3*(1-pi(i)))/((1-w)*l*(p1+(p2+p3)*pi(i))*gG+w*(1-l)*p3*(1-pi(i))*gB);
    R(3,2,i)=r*((1-w)*l*(p1+(p2+p3)*pi(i))+w*(1-l)*p3*(1-pi(i)))/((1-w)*l*(p1+(p2+p3)*pi(i))*gG+w*(1-l)*p3*(1-pi(i))*gB);
end

%pi cutoffs
phat=zeros(3,2);
phat(1,1)=-(w*l*p1*(y*gG-D*r)+(1-w)*(1-l)*p3*(y*gB-D*r))/(w*l*(p2+p3)*(y*gG-D*r)+(1-w)*(1-l)*(-p3)*(y*gB-D*r));
phat(3,1)=-(w*l*p3*(y*gG-D*r)+(1-w)*(1-l)*p1*(y*gB-D*r))/(w*l*(-p3)*(y*gG-D*r)+(1-w)*(1-l)*(p2+p3)*(y*gB-D*r));
phat(1,2)=-((1-w)*l*p1*(y*gG-D*r)+w*(1-l)*p3*(y*gB-D*r))/((1-w)*l*(p2+p3)*(y*gG-D*r)+w*(1-l)*(-p3)*(y*gB-D*r));
phat(3,2)=-((1-w)*l*p3*(y*gG-D*r)+w*(1-l)*p1*(y*gB-D*r))/((1-w)*l*(-p3)*(y*gG-D*r)+w*(1-l)*(p2+p3)*(y*gB-D*r));


%individual terms of marginal benefit function
Hx=zeros(3,lp);
Lx=zeros(3,lp);
for i=1:lp
    Hx(1,i)=((w*l*(p2+p3)*gG+(1-w)*(1-l)*(-p3)*gB)/(w*l+(1-w)*(1-l)))*(y-D*R(1,1,i));
    Hx(2,i)=((w*l*(-p2)*gG+(1-w)*(1-l)*(-p2)*gB)/(w*l+(1-w)*(1-l)))*(y-D*R(2,1,1));
    Hx(3,i)=((w*l*(-p3)*gG+(1-w)*(1-l)*(p2+p3)*gB)/(w*l+(1-w)*(1-l)))*(y-D*R(3,1,i));
    Lx(1,i)=(((1-w)*l*(p2+p3)*gG+w*(1-l)*(-p3)*gB)/((1-w)*l+w*(1-l)))*(y-D*R(1,2,i));
    Lx(2,i)=(((1-w)*l*(-p2)*gG+w*(1-l)*(-p2)*gB)/((1-w)*l+w*(1-l)))*(y-D*R(2,2,1));
    Lx(3,i)=(((1-w)*l*(-p3)*gG+w*(1-l)*(p2+p3)*gB)/((1-w)*l+w*(1-l)))*(y-D*R(3,2,i));
end

%check conditions to construct marginal benefit function
%omega > cutoff for H, omega < cutoff for L
%pi > cutoff for A, pi < cutoff for C
xH=zeros(1,lp);
xL=zeros(1,lp);
for i=1:lp
    %High
    if w>wbar(1,1,i) && pi(i)>phat(1,1),    ...A?
        xH(i)=Hx(1,i);  end
    if w>wbar(2,1,i),                       ...B?
        xH(i)=xH(i)+Hx(2,i); end
    if w>wbar(3,1,i) && pi(i)<phat(3,1),    ...C?
        xH(i)=xH(i)+Hx(3,i); end
    %Low
    if w<wbar(1,2,i) && pi(i)>phat(1,2),    ...A?
        xL(i)=Lx(1,i);  end
    if w<wbar(2,1,i),                       ...B?
        xL(i)=xL(i)+Lx(2,i); end
    if w<wbar(3,2,i) && pi(i)<phat(3,2),    ...C?
        xL(i)=xL(i)+Lx(3,i); end
end
    
%firm's optimal pi
piH=zeros(1,lp);
piL=zeros(1,lp);
for i=1:lp
    piH(i)=1-xH(i)^(-0.5);
    piL(i)=1-xL(i)^(-0.5);
end

disp(['pi cutoff for A,H:    ' num2str(phat(1,1))])
disp(['pi cutoff for C,H:    ' num2str(phat(3,1))])
disp(['pi cutoff for A,L:    ' num2str(phat(1,2))])
disp(['pi cutoff for C,L:    ' num2str(phat(3,2))])
disp(['omega cutoff for B,H: ' num2str(wbar(2,1,1))])
disp(['omega cutoff for B,L: ' num2str(wbar(2,2,1))])

figure(1)
plot(pi,pi,pi,xH)
if phat(1,1)>0, line([phat(1,1) phat(1,1)],[0 1]), end
if phat(3,1)>0, line([phat(3,1) phat(3,1)],[0 1]), end
set(gca,'YLim',[0 1])

figure(2)
plot(pi,pi,pi,xL)
if phat(1,2)>0, line([phat(1,2) phat(1,2)],[0 1]), end
if phat(3,2)>0, line([phat(3,2) phat(3,2)],[0 1]), end
set(gca,'YLim',[0 1])
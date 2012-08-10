tic

gamma_G=0.90;    ...prob. pays off if good
gamma_B=0.50;    ...prob. pays off if bad
y=5;             ...output
D=1;             ...investment size
r=1.01;          ...risk free rate
p1=0.5;          ...parameters for rating prod. function
p2=0.3;          ...see above
p3=0.2;

par=[gamma_G;   %1
     gamma_B;   %2
     y;         %3
     D;         %4
     r;         %5
     p1;        %6
     p2;        %7
     p3];       %8

w=0.5:0.001:1;
lw=length(w);

piH=zeros(1,lw);
piL=zeros(1,lw);
spreadA=zeros(1,lw);
spreadB=zeros(1,lw);
spreadC=zeros(1,lw);
MBh=zeros(1,lw);
MBl=zeros(1,lw);
MCh=zeros(1,lw);
MCl=zeros(1,lw);
probAgivH=zeros(1,lw);
probBgivH=zeros(1,lw);
probCgivH=zeros(1,lw);
probAgivL=zeros(1,lw);
probBgivL=zeros(1,lw);
probCgivL=zeros(1,lw);

for k=1:lw
    [piH(k),piL(k),spreadA(k),spreadB(k),spreadC(k),MBh(k),MBl(k),MCh(k),MCl(k)]=spreadcalc(w(k),par);
end

probAgivH=par(6)+(par(7)+par(8))*piH;
probBgivH=par(7)*(1-piH);
probCgivH=par(8)*(1-piH);
probAgivL=par(8)*(1-piL);
probBgivL=par(7)*(1-piL);
probCgivL=par(6)+(par(7)+par(8))*piL;

figure(1)
subplot(2,2,1,'replace')
plot(w,piH,'--',w,piL)
title('Optimal investment in ratings')
xlabel('\omega')
ylabel('\pi^*_\nu')
legend('\pi^*_H','\pi^*_L','Location','best')

subplot(2,2,2,'replace')
plot(w,spreadA,'--g',w,spreadB,'-b',w,spreadC,'r')
title('Interest rate spreads (R(h,L)-R(h,H))')
xlabel('\omega')
legend('A','B','C','Location','best')

subplot(2,2,3,'replace')
plot(w,MBh,'--',w,MBl)
title('Marginal Benefit')
xlabel('\omega')
legend('MB H','MB L','Location','best')

subplot(2,2,4,'replace')
plot(w,MCh,'--',w,MCl)
title('Marginal Cost')
xlabel('\omega')
ylim=([-0.5 1]);
legend('MC H','MC L','Location','best')

figure(2)
hold on
h1=plot(w,probAgivH,'b');
h1b=plot(w(1:10:end),probAgivH(1:10:end),'xb');

h2=plot(w,probBgivH,'c');
h2b=plot(w(1:10:end),probBgivH(1:10:end),'xc');

h3=plot(w,probCgivH,'m');
h3b=plot(w(1:10:end),probCgivH(1:10:end),'mx');

h4=plot(w,probCgivL,'r');
h4b=plot(w(1:10:end),probCgivL(1:10:end),'xr');

h5=plot(w,probBgivL,'g');
h5b=plot(w(1:10:end),probBgivL(1:10:end),'gx');

h6=plot(w,probAgivL,'y');
h6b=plot(w(1:10:end),probAgivL(1:10:end),'yx');

legend([h1 h2 h3 h4 h5 h6],'AH','BH','CH','CL','BL','AL','Location','best')
hold off

toc
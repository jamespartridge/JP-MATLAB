tic

w=0.5:0.01:1;
lw=length(w);

piH=zeros(1,lw);
piL=zeros(1,lw);
spreadA=zeros(1,lw);
spreadB=zeros(1,lw);
spreadC=zeros(1,lw);
MBh=zeros(1,lw);
MBl=zeros(1,lw);


for k=1:lw
    [piH(k),piL(k),spreadA(k),spreadB(k),spreadC(k),MBh(k),MBl(k)]=spreadcalc(w(k));
    k
end

figure(1)
plot(w,piH,'g',w,piL,'b')
legend('pi*H','pi*L','Location','best')

figure(2)
plot(w,spreadA,'g',w,spreadB,'b',w,spreadC,'r')
legend('A','B','C','Location','best')

figure(3)
plot(w,MBh,w,MBl)
legend('MB H','MB L','Location','best')

toc
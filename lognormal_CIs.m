clear all;
tic
alpha=0.05;

N=zeros(4,5);
S=zeros(4,5);
lower=zeros(4,5);
cov=zeros(4,5);
upper=zeros(4,5);

% Matrix S is the collected standard deviations where a row is a rating;
% class (All, HI, LO, SG) and a column is a period (90, 94, 98, 02, 06);
% Matrix N is the collected observations with the same ordering as S;
S=  [0.5781	0.8247	0.7725	0.8014	0.6988;
     0.4178	0.3977	0.6272	0.5421	0.6532;
     0.4366	0.4381	0.5139	0.4993	0.5633;
     0.5284	0.4509	0.3866	0.3906	0.4681];

N=  [2797	4489	5168	4581	2954;
     424	433     531     553     298;
     1666	2608	2291	1490	1389;
     678	1425	2329	2524	1263];

% Starting guess;
x0=[0,1];
for i=1:5
    for j=1:4
        n=N(j,i);
        s=S(j,i);
        [x,fval] = runconfint(n,s,alpha,x0);
        lower(j,i)=sqrt(exp(x(1))-1);
        cov(j,i)=sqrt(exp(s)-1);
        upper(j,i)=sqrt(exp(x(2))-1);
    end
end
toc

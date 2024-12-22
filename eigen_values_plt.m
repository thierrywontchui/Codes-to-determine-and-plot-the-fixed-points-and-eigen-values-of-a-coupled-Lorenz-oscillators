% This code determine and plot the maximum eigen values of the Jacobian
% matrix of the coupled Lorenz system with coupling strength c=Ex=Ez.

clear all
close all 
clc

alpha=10; beta=24.76; gama=8/3;

syms x1 y1 z1 x2 y2 z2 c 

eqn1=10*(y1-x1)+c*(x2-x1)==0;
eqn2=24.76*x1-y1-x1*z1==0;
eqn3=-(8/3)*z1+x1*y1+c*(z2-z1)==0;
eqn4=10*(y2-x2)+c*(x1-x2)==0;
eqn5=24.76*x2-y2-x2*z2==0;
eqn6=-(8/3)*z2+x2*y2+c*(z1-z2)==0;

S=vpasolve([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6],[x1,y1,z1,x2,y2,z2]);
S.x1;
S.y1;
S.z1;
S.x2;
S.y2;
S.z2;
solu=[S.x1,S.y1,S.z1,S.x2,S.y2,S.z2];

c=0.0000001:((12.804114685453-0.0000001)/40)+0:12.804114685453; % first semgment of the coupling strength range. Plotting in piece enable us to change the solution branch colors depending on their nature.
% c=0.0000001:(0.00444-0.0000001)/20:0.00444;
% c=12.804114685453:((20-12.804114685453)/22)+0:20;             % second segment of the coupling strength range.
% cp=c;
lamda1=zeros(length(c),1);
lamda2=zeros(length(c),1);
lamda3=zeros(length(c),1);
lamda4=zeros(length(c),1);
lamda5=zeros(length(c),1);
lamda6=zeros(length(c),1);
lamda7=zeros(length(c),1);
lamda8=zeros(length(c),1);
lamda9=zeros(length(c),1);

for k=1:length(c)
    sol=subs(solu,c(k));
    lambd=zeros(9,1);
    for i=1:9
        vsol=sol(i,:);
        x1=vsol(1,1);
        y1=vsol(1,2);
        z1=vsol(1,3);       
        x2=vsol(1,4);
        y2=vsol(1,5);
        z2=vsol(1,6);
        J=[-alpha-c(k) alpha 0 c(k) 0 0;
            beta-z1 -1 -x1 0 0 0;
            y1 x1 -gama-c(k) 0 0 c(k);
            c(k) 0 0 -alpha-c(k) alpha 0;
            0 0 0 beta-z2 -1 -x2;
            0 0 c(k) y2 x2 -gama-c(k)];
        lambda=eig(J);
        lambdap=max(real(lambda));
        lambd(i,1)=lambdap;
    end
    lamda1(k,1)=lambd(1,1);
    lamda2(k,1)=lambd(2,1);
    lamda3(k,1)=lambd(3,1);
    lamda4(k,1)=lambd(4,1);
    lamda5(k,1)=lambd(5,1);
    lamda6(k,1)=lambd(6,1);
    lamda7(k,1)=lambd(7,1);
    lamda8(k,1)=lambd(8,1);
    lamda9(k,1)=lambd(9,1);
end

% save EV cp lamda1p lamda2p lamda3p lamda4p lamda5p lamda6p lamda7p lamda8p lamda9p   %% These save maximum eigen value solution of J when c=12.8--20. 
load EV                                              % EV contains the eigen value solutions of J for c=12.8--20.
plot(c,lamda1, '-.g', 'LineWidth',1.4)
hold on
plot(c,lamda2, '-.g', 'LineWidth',1.4)
plot(c,lamda3, '-b', 'LineWidth',1.5)
plot(c,lamda4, '-b', 'LineWidth',1.5)
plot(c,lamda5, '.r', 'MarkerSize',8.0)
plot(c,lamda6, '--r', 'LineWidth',1.5)
plot(c,lamda7, '--r', 'LineWidth',1.5)
plot(c,lamda8, '--r', 'LineWidth',1.5)
plot(c,lamda9, '--r', 'LineWidth',1.5)

plot(cp,lamda1p, '-.g', 'LineWidth',1.4)
plot(cp,lamda2p, '-.g', 'LineWidth',1.4)
plot(cp,lamda3p, '-r', 'LineWidth',1.5)
plot(cp,lamda4p, '-r', 'LineWidth',1.5)
plot(cp,lamda5p, '.r', 'MarkerSize',8.0)

set(gca,'FontSize',20,'FontName','Times','fontangle','italic')
xlabel([char(949)])
ylabel('max(Re(\lambda))_ ^ ')
set(gca,'YTick',-5:5:15)
set(gca,'XTick',0:5:20)
text(18.5,14.0,'(c)','FontSize',16)
xlim([0 20])
ylim([-5 15])
grid on
hold off        
    
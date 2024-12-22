% This code  determine the symbolic fixed point solutions of the coupled
% Lorenz system when the interaction is simultaneously done through both the -x and -z
% directions.

clear all
close all 
clc

syms x1 y1 z1 x2 y2 z2 Ex Ez 

eqn1=10*(y1-x1)+Ex*(x2-x1)==0;
eqn2=24.76*x1-y1-x1*z1==0;
eqn3=-(8/3)*z1+x1*y1+Ez*(z2-z1)==0;
eqn4=10*(y2-x2)+Ex*(x1-x2)==0;
eqn5=24.76*x2-y2-x2*z2==0;
eqn6=-(8/3)*z2+x2*y2+Ez*(z1-z2)==0;

% S=vpasolve([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6],[x1,y1,z1,x2,y2,z2]); %% This solver can be used (uncomment) if the fixed points solutions have to be output on decimal form. However the row solutions might change position in the matrix solutions. 

S=solve([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6],[x1,y1,z1,x2,y2,z2]);

S.x1;
S.y1;
S.z1;
S.x2;
S.y2;
S.z2;

% the output "solutions" is the matrix solutions 9*6  of the system. The 9 rows
% represent the nine possible solutions of the coupled system and the 6
% columns represent solution x1, y1, z1, x2, y2, z2, respectively.

solutions=[S.x1,S.y1,S.z1,S.x2,S.y2,S.z2]

S1= solutions(3,:)
S2= solutions(2,:)
S3= solutions(1,:)
S4= solutions(4,:)
S5= solutions(5,:)
S6= solutions(6,:)
S7= solutions(8,:)
S8= solutions(7,:)
S9= solutions(9,:)


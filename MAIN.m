%% NOTICE:
% S: BASIC VARIABLES,[rho u p]
% Q: CONSERVATIVE VARIABLES,[rho rho*u rho*E]
% F: FLUX VARIABLES,[rho*u u^2+p rho*u*H]
clc;
clear;
path(path,'.\tools');
path(path,'.\schemes');

%% fluid properties
global gamma R;
gamma=1.4;
R=287;

%% common settings
N=100;% grid number
CFL=0.8;% CFL number
x=linspace(0,10,N)';% grid
delta_x=(x(end)-x(1))/(N-1);% grid spacing
convergence=1e-6;

%% boundary conditions
M_infty=1.5;
p_infty=47892.4;
rho_infty=1.2218;
u_infty=M_infty*sqrt(gamma*p_infty/rho_infty);
S_infty=[rho_infty,u_infty,p_infty];
Q_infty=S2Q(S_infty);

u_out=119;% if it's a subsonic boundary

%% initialization
S=zeros(N,3);% initial field
S(:,1)=rho_infty;
S(:,2)=u_infty;
S(:,3)=p_infty;

%% solving
supersonicSolution=numericalSolver(x,S,CFL,convergence,Q_infty);
subsonicSolution=numericalSolver(x,S,CFL,convergence,Q_infty,u_out);

%% ploting
ploting;







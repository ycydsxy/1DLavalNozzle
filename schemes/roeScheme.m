function F = roeScheme(Q_l,Q_r)
% Computing F_i+1/2 Using Roe Scheme
%   Q_l: [rho rho*u rho*E] on the left side
%   Q_r: [rho rho*u rho*E] on the right side
%   F: [rho*u u^2+p rho*u*H] 

global gamma;

S_l=Q2S(Q_l);
S_r=Q2S(Q_r);
F_l=S2F(S_l);
F_r=S2F(S_r);

rho_1=S_l(1);
u_1=S_l(2);
p_1=S_l(3);
h_1=gamma/(gamma-1)*p_1/rho_1+0.5*u_1^2;
rho_2=S_r(1);
u_2=S_r(2);
p_2=S_r(3);
h_2=gamma/(gamma-1)*p_2/rho_2+0.5*u_2^2;

D=sqrt(rho_2/rho_1);
u=(D*u_2+u_1)/(D+1);
h=(D*h_2+h_1)/(D+1);
c=sqrt((gamma-1)*(h-0.5*u^2));

Left=[1,1,1;u,u-c,u+c;0.5*u^2,h-c*u,h+c*u];
Gamma=abs(diag([u,u-c,u+c]));
b_t1=(gamma-1)/c^2*(0.5*u^2*(Q_r(1)-Q_l(1))-u*(Q_r(2)-Q_l(2))+(Q_r(3)-Q_l(3)));
b_t2=1/c*(u*(Q_r(1)-Q_l(1))-(Q_r(2)-Q_l(2)));
b_t3=(Q_r(1)-Q_l(1))-b_t1;
Right=[b_t3;0.5*(b_t1+b_t2);0.5*(b_t1-b_t2)];

F_d=-0.5*Left*Gamma*Right;

F=0.5*(F_l+F_r)+F_d';
end


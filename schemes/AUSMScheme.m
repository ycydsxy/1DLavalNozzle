function F = AUSMScheme(Q_l,Q_r)
% Computing F_i+1/2 Using AUSM Scheme
%   Q_l: [rho rho*u rho*E] on the left side
%   Q_r: [rho rho*u rho*E] on the right side
%   F: [rho*u u^2+p rho*u*H] 

global gamma;

S_l=Q2S(Q_l);
S_r=Q2S(Q_r);

rho_l=S_l(1);
u_l=S_l(2);
p_l=S_l(3);
h_l=gamma/(gamma-1)*p_l/rho_l+0.5*u_l^2;
c_l=sqrt(gamma*p_l/rho_l);
Ma_l=u_l/c_l;

rho_r=S_r(1);
u_r=S_r(2);
p_r=S_r(3);
h_r=gamma/(gamma-1)*p_r/rho_r+0.5*u_r^2;
c_r=sqrt(gamma*p_r/rho_r);
Ma_r=u_r/c_r;

if abs(Ma_l)<=1
    lambda_l=0.25*(Ma_l+1)^2;
else
    lambda_l=0.5*(Ma_l+abs(Ma_l));
end

if abs(Ma_r)<=1
    lambda_r=-0.25*(Ma_r-1)^2;
else
    lambda_r=0.5*(Ma_r-abs(Ma_r));
end

Ma_c=lambda_l+lambda_r;

if Ma_c>=0
    c_c=c_l;
else
    c_c=c_r;
end

m_c=c_c*Ma_c;

if m_c>=0
   F_c=m_c*[rho_l,rho_l*u_l,rho_l*h_l]; 
else
   F_c=m_c*[rho_r,rho_r*u_r,rho_r*h_r]; 
end

if abs(Ma_l)>1
   psi_l=0.5*(Ma_l+abs(Ma_l))/Ma_l;
else
   psi_l=0.25*(Ma_l+1)^2*(2-Ma_l);
end

if abs(Ma_r)>1
   psi_r=0.5*(Ma_r-abs(Ma_r))/Ma_r;
else
   psi_r=0.25*(Ma_r-1)^2*(2+Ma_r);
end

P_c=psi_l*[0,p_l,0]+psi_r*[0,p_r,0];

F=F_c+P_c;
end


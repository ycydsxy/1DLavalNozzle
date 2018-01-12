function [ F ] = S2F( S )
% transform S to F
%   S: [rho u p]
%   F: [rho*u rho*u^2+p rho*u*H]

global gamma;

rho=S(:,1);
u=S(:,2);
p=S(:,3);

F=[rho.*u,rho.*u.^2+p,p.*u*gamma/(gamma-1)+0.5*rho.*u.^3];

end


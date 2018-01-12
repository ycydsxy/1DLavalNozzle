function [ Q ] = S2Q( S )
% transform S to Q
%   S: [rho u p]
%   Q: [rho rho*u rho*E]

global gamma;

rho=S(:,1);
u=S(:,2);
p=S(:,3);

Q=[rho,rho.*u,p/(gamma-1)+0.5*rho.*u.^2];

end


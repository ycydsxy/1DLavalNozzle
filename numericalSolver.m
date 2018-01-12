function result = numericalSolver(x,S,CFL,convergence,Q_in,u_out)
% Solver Using Different Schemes
%   x: value of coordinate x
%   S: [rho u p]
%   CFL: CFL number
%   convergence: convergence criteria
%   Q_in: inlet boundary conditions
%   u_out[optional]: outlet boundary conditions(if subsonic)
%   result: [rho u p]
tic;

global gamma;

if nargin<6
    boundaryType=1;% supersonic boundary conditions
    disp('supersonic boundary conditions.')
else
    boundaryType=2;% subsonic boundary conditions
    disp('subsonic boundary conditions.')
end

N=length(x);
delta_x=(x(end)-x(1))/(N-1);
Q=S2Q(S);
F=zeros(N-1,3);
currentStep=0;
maxSteps=1e5;
history=zeros(maxSteps,1);

while(1)
    currentStep=currentStep+1;
    S_last=S;
    delta_t=CFL*delta_x/max(abs(S(:,2))+sqrt(gamma*S(:,3)./S(:,1)));% timestep
    
    for i=1:N-1 % computing F_i+1/2
        F(i,:)=roeScheme(Q(i,:),Q(i+1,:));
    end
    
    source=[zeros(N,1),S(:,3),zeros(N,1)]-S2F(S);% computing source
    source_temp=0.347*0.8.*(sech(0.8*x-4)).^2./(1.398+0.347*tanh(0.8*x-4));
    source(:,1)=source(:,1).*source_temp;
    source(:,2)=source(:,2).*source_temp;
    source(:,3)=source(:,3).*source_temp;
    
    for i=2:N-1 % computing Q_i^(n+1) using FISRT ORDER time scheme
        Q(i,:)=Q(i,:)-delta_t/delta_x*(F(i,:)-F(i-1,:))+delta_t*source(i,:);
    end
    
    switch boundaryType
        case 1% supersonic
            Q(1,:)=Q_in;% inlet boundary, needn't outlet boundary
            Q(N,:)=Q(N-1,:);
        case 2% subsonic
            Q(1,:)=Q_in;% inlet boundary
            S_out=Q2S(Q(N-1,:));
            S_out(2)=u_out;
            Q(N,:)=S2Q(S_out);
        otherwise
            error('The Boundary You Choose ISNOT Supported!')
    end
    
    S=Q2S(Q);
    
    residual=max(max(abs(S-S_last)));
    history(currentStep)=residual;
    if currentStep>maxSteps % stoping criteria
        disp('WARNING: maxSteps reached!');
        break;
    end
    if residual<convergence % stoping criteria
        break;
    end
end

toc

figure;
plot(history(history~=0),'Color','r');
set(gca,'YScale','log');

result=S;

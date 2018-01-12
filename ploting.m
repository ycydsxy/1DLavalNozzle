%% ploting
ylabels={'\rho','u','p'};

figure
for j=1:3
    subplot(3,1,j);
    if j==1
        title(sprintf('SUPERSONIC CFL=%.2f,N=%d',CFL,N));
    end
    hold on
    plot(x,supersonicSolution(:,j),'o');
    
    xlabel('x');
    ylabel(ylabels{j});
end

figure
for j=1:3
    subplot(3,1,j);
    if j==1
        title(sprintf('SUBSONIC CFL=%.2f,N=%d',CFL,N));
    end
    hold on
    plot(x,subsonicSolution(:,j),'o');
    
    xlabel('x');
    ylabel(ylabels{j});
end
function lorenz_diagram()
sigma = 10;
beta = 8/3;
% the function giving the asymptotically stable equilibrium C_{+}
cplus = @(rho)[sqrt(beta*(rho-1)),sqrt(beta*(rho-1)),rho-1]';

% assign the desired range of rho and step
rmin = 1.05; % must be > 1
rmax = 30;% 350;
rstep = 0.1;

% set up figure
fig = 3;
figure(fig); clf;
grid;
hold on;
set(gca,'FontSize',20);
xlabel('\rho','FontSize',20);
ylabel('x_1','FontSize',20);
Mksize = 6; % markersize
col = [1,  0.6,0.6; % pink
       0.6,0.6,0.6; % grey
       1,  0,  0;   % red
       0,  0,  0];  % black 
k = 0;
for rho = rmin : rstep : rmax
    lo = @(t,x) [-sigma*x(1) + sigma*x(2); x(1).*(rho - x(3)) - x(2); -beta*x(3) + x(1).*x(2)];
    Cp = cplus(rho);
    sq = sqrt(beta*(rho-1));
    myevents = @(t,y)crossplane(y,Cp); % horizontal plane z = rho - 1
    options = odeset('AbsTol',1e-6,'RelTol',1e-6);
    % find xi = the unstable direction for the origin
    J = [-sigma,sigma,0;rho,-1,0;0,0,-beta];
    [V,E] = eig(J);
    [esort,isort] = sort(diag(E),'descend');
    V = V(:,isort);
    if V(1,1) < 0
        V(:,1) = -V(:,1);
    end
    y0 = 0.01*V(:,1); % trajectories emanate from the origin along the unstable direction
    options = odeset('AbsTol',1e-6,'RelTol',1e-6,'Events',myevents);
    % record all crosses of the horizontal plane z = rho - 1
    [T,Y,Te,Ye,Ie] = ode45(lo,[0,200],y0,options);
    k = k + 1;
    if ~isempty(Te)
        figure(fig);
        plot(rho*ones(size(Ye,1),1),Ye(:,1),'.','Markersize',Mksize,'color',col(1,:));
        plot(rho*ones(size(Ye,1),1),-Ye(:,1),'.','Markersize',Mksize,'color',col(2,:));
        if mod(k,100) == 0
            drawnow;
        end
    else
        plot(rho,Cp(1),'.','Markersize',Mksize,'color',col(1,:));
        plot(rho,-Cp(1),'.','Markersize',Mksize,'color',col(2,:));
        if mod(k,100) == 0
            drawnow;
        end
    end 
    % continue
    % at this point, the characteristic settles on the attractor
    y0 = Y(end,:);
    [T,Y,Te,Ye,Ie] = ode45(lo,[0,200],y0,options);
    % record all crosses of the horizontal plane z = rho - 1
    k = k + 1;
    if ~isempty(Te)
        figure(fig);
        plot(rho*ones(size(Ye,1),1),Ye(:,1),'.','Markersize',Mksize,'color',col(3,:));
        plot(rho*ones(size(Ye,1),1),-Ye(:,1),'.','Markersize',Mksize,'color',col(4,:));
        if mod(k,100) == 0
            drawnow;
        end
    else
        if norm(y0' - Cp) < norm(y0' -[-Cp(1:2);Cp(3)])
            ss = 1;
        else
            ss = -1;
        end
        plot(rho,ss*Cp(1),'.','Markersize',Mksize,'color',col(3,:));
        plot(rho,-ss*Cp(1),'.','Markersize',Mksize,'color',col(4,:));
        if mod(k,100) == 0
            drawnow;
        end
    end    

end
% plot vertical lines corresponding to the critical values of rho
ymax = 24;
plot([13.926,13.926],[-ymax,ymax],'--','Linewidth',2,'color',[0,0.6,0])
plot([24.06,24.06],[-ymax,ymax],'--','Linewidth',2,'color',[0,0.6,0])
plot([24.74,24.74],[-ymax,ymax],'--','Linewidth',2,'color',[0,0.6,0])

end

%%

function [position,isterminal,direction] = crossplane(y,ye)
position = ye(3) - y(3); % The value that we want to be zero
isterminal = 0;  % Do not halt integration 
direction = 0;   
end



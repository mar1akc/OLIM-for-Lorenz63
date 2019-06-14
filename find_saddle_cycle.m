function [Y1loop,l1loop,l1] = find_saddle_cycle(rho)
% Finds saddle limit cycle for 13.926 < rho < 24.74 
global lo b V x0 ye ye1
sigma = 10;
beta = 8/3;
% Pick figure number at which you want to plot
% the equilibria C+ and C-, the cycles gamma+ and gamma-, and the mesh

%% First find the saddle cycle

ye = [sqrt(beta*(rho-1));sqrt(beta*(rho-1));rho-1]; % equilibrium
ye1 = [-sqrt(beta*(rho-1));-sqrt(beta*(rho-1));rho-1]; % equilibrium
A = [-sigma,sigma,0;1,-1,-ye(1);ye(1),ye(1),-beta];
U = inv(sylvester(A,A',-2*eye(3)));
[eV,eE] = eig(U);
eval = diag(eE);
[eval,isort] = sort(eval,'descend');
evec = eV(:,isort); % the largest eigenpair goes first
v1 = evec(:,1);

% find the initial guess for finding limit cycle
lo = @(t,a) [-sigma*a(1) + sigma*a(2); a(1)*(rho - a(3)) - a(2); -beta*a(3) + a(1)*a(2)];
options1 = odeset('AbsTol',1e-12,'RelTol',1e-12);
options2 = odeset('AbsTol',1e-12,'RelTol',1e-12,'Events',@myevents);

flag = 0;
step = 0.1;
j = 1;
% plot3(ye(1),ye(2),ye(3),'r.','Markersize',30);
% plot3(-ye(1),-ye(2),ye(3),'r.','Markersize',30);
% view(3);
% drawnow;

while flag == 0
    x0 =  ye*step*j;  %ye + (r + j*step)*evec(:,2);
    [~,Y] = ode45(lo,[0,20],x0,options1);
    if min(Y(:,1)) > 0
        flag = 1;
    end
    j = j + 1;
end    
b = lo(0,x0); b = b/norm(b);

A = [-sigma,sigma,0;1,-1,-ye(1);ye(1),ye(1),-beta];
U = inv(sylvester(A,A',-2*eye(3)));
[eV,eE] = eig(U);
eval = diag(eE);
[eval,isort] = sort(eval,'descend');
evec = eV(:,isort); % the largest eigenpair goes first
v1 = evec(:,1);

aux = b + v1;

v1 = cross(b,aux); v1 = v1/norm(v1);
v2 = cross(v1,b); v2 = v2/norm(v2);
% {v1, v2} = 0NB in the plane normal to b and passing through x0
V = [v1, v2];

x = x0; 

iter = 0;
maxiter = 1000;
pp = V'*poincare(x)';
while norm(pp) > 1e-12 && iter <= maxiter
    J = jacobian(x);
    p=-V*(J\pp);
    if norm(p) > 1
        p = p/norm(p);
    end
    x = x + p;
    iter=iter+1;
    fprintf('i = %d, norm = %d, x = [%d,%d,%d]\n',iter,norm(poincare(x)),x(1),x(2),x(3));
    pp = V'*poincare(x)';
end

[~,X] = ode45(lo,[0 1e-3],x,options1);
[~,Y2,~,~,~] = ode45(lo,[0 6],X(end,:),options2);
Y2 = [X;Y2];
% plot3(Y2(:,1),Y2(:,2),Y2(:,3),'Linewidth',4,'color','r');
% plot3(-Y2(:,1),-Y2(:,2),Y2(:,3),'Linewidth',4,'color','r');
% drawnow;
% view(3)
fname = sprintf('LorenzLimitCycle_rho_%.2f.mat',rho);
save(fname,'Y2');
Y1loop = Y2;

[Y1loop,l1loop,l1] = arclength(Y1loop);
fprintf('length(Y1loop) = %d\n',l1);
Nloop = size(Y1loop,1);
fprintf('Nloop = %d\n',Nloop - 1);
options3 = odeset('AbsTol',1e-12,'RelTol',1e-12,'Events',@myevents3);

% 
% daspect([1,1,1])
% set(gca,'FontSize',20);
% xlabel('x_1');
% ylabel('x_2');
% zlabel('x_3');
% 

end

%%
function dd = poincare(x)
global lo
options1 = odeset('AbsTol',1e-12,'RelTol',1e-12);
options2 = odeset('AbsTol',1e-12,'RelTol',1e-12,'Events',@myevents);
[~,X] = ode45(lo,[0 1e-3],x,options1);
[~,Y,~,~,~] = ode45(lo,[0 6],X(end,:),options2);
Y = [X;Y];
dd = Y(end,:)-x';
end

function g = jacobian(x)
global V
h = 0.01;
g = zeros(2);
g11 = poincare(x + V(:,1)*h);
g12 = poincare(x - V(:,1)*h);
g21 = poincare(x + V(:,2)*h);
g22 = poincare(x - V(:,2)*h);
g(:,1) = (g11 - g12)*V/(2*h);
g(:,2) = (g21 - g22)*V/(2*h);
end

function output=func(x)
output=0.5*norm(poincare(x))^2;     
end

%%

function [position,isterminal,direction] = myevents(~,y)
global  b x0
position = (y - x0)'*b; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 1;   
end

%%

function [position,isterminal,direction] = myevents3(~,y)
global  ye ye1
position = min([norm(y - ye),norm(y - ye1)]) - 1e-2; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 0;   
end

%%
function [Y0,lxyz,lY] = arclength(Y0)
% rows of Y0 must have three entries: x, y ,z
% Y0 is cleaned from repeated points
% lxyz = the arclength parametrization of Y0
% lY = length of Y0
dx = Y0 - circshift(Y0,[1,0]);
dx(1,:) = [0,0,0];
dl = sqrt(sum(dx.^2,2));
lxyz = cumsum(dl);
lY = lxyz(end);
[lxyz,iu,~] = unique(lxyz);
Yaux = Y0(iu,:);
clear Y0
Y0 = Yaux;
end


%%
function [position,isterminal,direction] = MyCross(y,a,y0)
position = (y - y0)'*a; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 1;   
end



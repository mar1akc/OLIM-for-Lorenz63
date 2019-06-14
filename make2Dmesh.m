function make2Dmesh()
% Finds limit cycle for 13.926 < rho < 24.74 and generates a 2D radial mesh
% between the stable equilibrium and the limit cycle
global lo b V x0
rho = 20.0;
sigma = 10;
beta = 8/3;
% The size of the radial mesh will be Nstart by Ncut
% If you want to plot the mesh, mare Nstart = 40, and Ncut = 72
Nstart = 40; %2001;
Ncut = 72; %7201;
% If plot_mesh_flag == 1, the mesh will be visualised,
plot_mesh_flag = 1;
% Pick figure number at which you want to plot
% the equilibria C+ and C-, the cycles gamma+ and gamma-, and the mesh
fig = 6;
%% File names where the mesh will be saved

% The mat file with the mesh and all parameters
fname_mat = sprintf('DemoMesh2Drho%.2f.mat',rho);

% The txt file with mesh size data
fname1 = sprintf('DemoMesh2Ddata_rho%.2f.txt',rho);

% The txt file with mesh points
fname2 = sprintf('DemoMesh2Drho%.2f.txt',rho);

%% First find the saddle cycle
figure(fig); clf; hold on;grid

ye = [sqrt(beta*(rho-1));sqrt(beta*(rho-1));rho-1]; % equilibrium
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
plot3(ye(1),ye(2),ye(3),'r.','Markersize',30);
plot3(-ye(1),-ye(2),ye(3),'r.','Markersize',30);
view(3);
drawnow;
while flag == 0
    x0 =  ye*step*j;  %ye + (r + j*step)*evec(:,2);
    [~,Y] = ode45(lo,[0,20],x0,options1);
    if min(Y(:,1)) > 0
        flag = 1;
    end
    j = j + 1;
end    
b = lo(0,x0); b = b/norm(b);
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
    iter = iter+1;
    fprintf('i = %d, norm = %d, x = [%d,%d,%d]\n',iter,norm(poincare(x)),x(1),x(2),x(3));
    pp = V'*poincare(x)';
end

[~,X] = ode45(lo,[0 1e-3],x,options1);
[~,Y2,~,~,~] = ode45(lo,[0 6],X(end,:),options2);
Y2 = [X;Y2];
plot3(Y2(:,1),Y2(:,2),Y2(:,3),'Linewidth',4,'color','r');
plot3(-Y2(:,1),-Y2(:,2),Y2(:,3),'Linewidth',4,'color','r');
drawnow;
view(3)
fname = sprintf('LorenzLimitCycle_rho_%.2f.mat',rho);
save(fname,'Y2');
Y1loop = Y2;

%% make a mesh
xs = x + 1e-3*(ye - x);
y0(1,:) = x';
k = 1;
    a = xs - ye;
    a = a/norm(a);
    aa = cross(a,lo(0,xs));
    aa = aa/norm(aa);
    a = cross(aa,a);
    events = @(t,y)MyCross(y,a,x);
    options2 = odeset('AbsTol',1e-12,'RelTol',1e-12,'refine',10,'Events',events);

while norm(ye - xs) > 1e-3
    [~,X] = ode45(lo,[0 1e-3],xs,options1);
    [~,Y,~,~,~] = ode45(lo,[0 6],X(end,:),options2);
    Y = [X;Y];
    xs = Y(end,:)';
    k = k + 1;
    y0(k,:) = Y(end,:);
end
k = k + 1;
y0(k,:) = ye;
[y0,lxyz,ly0] = arclength(y0);
%parametrize y0 uniformly
ystart = interp1(lxyz,y0,linspace(0,ly0,Nstart)');
[xmesh,ymesh,zmesh] = MakeMesh(ye,Y1loop,Nstart,Ncut,ystart(end:-1:1,:));

if plot_mesh_flag == 1
    plotmesh(xmesh,ymesh,zmesh);
    plotmesh(-xmesh,-ymesh,zmesh)
end

daspect([1,1,1])
set(gca,'FontSize',20);
xlabel('x_1');
ylabel('x_2');
zlabel('x_3');
view(3);

%% record the mesh
save(fname_mat,'xmesh','ymesh','zmesh','rho','sigma','beta','Ncut','Nstart');

fid = fopen(fname1,'w');
fprintf(fid,'%d\n%d\n',Ncut,Nstart);
fclose(fid);

fid = fopen(fname2,'w');
for j = 1 : Nstart
    fprintf(fid,'%.14e\t',xmesh(j,1 : end - 1));
    fprintf(fid,'\n');
end
for j = 1 : Nstart
    fprintf(fid,'%.14e\t',ymesh(j,1 : end - 1));
    fprintf(fid,'\n');
end
for j = 1 : Nstart
    fprintf(fid,'%.14e\t',zmesh(j,1 : end - 1));
    fprintf(fid,'\n');
end
fclose(fid);
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
%%
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

%%
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
function [xmesh,ymesh,zmesh] = MakeMesh(ye,Y1loop,Nstart,Ncut,ystart)
global lo 
[Y1loop,l1loop,l1] = arclength(Y1loop);
fprintf('length(Y1loop) = %d\n',l1);

gstart = linspace(0,1,Nstart)';
gcut = linspace(0,1,Ncut)';

% mesh points at the first layer all collapse into the initial equilibrium
xmesh = zeros(Nstart,Ncut);
xmesh(1,:) = ye(1)*ones(1,size(xmesh,2));
ymesh = zeros(Nstart,Ncut);
ymesh(1,:) = ye(2)*ones(1,size(xmesh,2));
zmesh = zeros(Nstart,Ncut);
zmesh(1,:) = ye(3)*ones(1,size(xmesh,2));

% mesh points at the last layer all lie at the limit cycle
Y1cutarray = interp1(l1loop,Y1loop,gcut*l1)';
xmesh(end,:) = Y1cutarray(1,:);
ymesh(end,:) = Y1cutarray(2,:);
zmesh(end,:) = Y1cutarray(3,:);

% mesh points of the first cut are found
xmesh(:,1) = ystart(:,1);
ymesh(:,1) = ystart(:,2);
zmesh(:,1) = ystart(:,3);

% mesh points of the last cut are the same as at the first one 
xmesh(:,Ncut) = ystart(:,1);
ymesh(:,Ncut) = ystart(:,2);
zmesh(:,Ncut) = ystart(:,3);

for k = 2 : Ncut - 1    
    fprintf('Cut # %i\n',k);
    Y1cut = Y1cutarray(:,k);
    bcut = lo(0,Y1cut);
    bcut = bcut/norm(bcut);
    a = Y1cut - ye;
    a = a/norm(a);
    aa = cross(a,bcut);
    aa = aa/norm(aa);
    b = cross(aa,a);
    events = @(t,y)MyCross(y,b,Y1cut);
    options1 = odeset('AbsTol',1e-12,'RelTol',1e-12,'refine',10,'Events',events);
    yend = zeros(Nstart,3);
    yend(1,:) = ye';
    yend(end,:) = Y1cut';
    for j = 2 : Nstart - 1
        y0 = ystart(j,:);
        [~,Y] = ode45(lo,[0,100],y0,options1);
        aux = Y(end,:);
        yend(j,:) = Y(end,:);
    end
    d = sqrt(sum((yend - ones(size(yend,1),1)*ye').^2,2));
    [dsort,isort] = sort(d,'ascend');
    yend = yend(isort,:);
    
    [yend,lyend,ly] = arclength(yend);
    ystart = interp1(lyend,yend,linspace(0,ly,Nstart)','spline');
    xmesh(:,k) = ystart(:,1);
    ymesh(:,k) = ystart(:,2);
    zmesh(:,k) = ystart(:,3);
end
end

%%
function [position,isterminal,direction] = MyCross(y,a,y0)
position = (y - y0)'*a; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 1;   
end

%%
function plotmesh(xmesh,ymesh,zmesh)
[Nstart,Ncut] = size(xmesh);
for k = 1 : Ncut
    plot3(xmesh(:,k),ymesh(:,k),zmesh(:,k));
end
for k = 1 : Nstart
    plot3(xmesh(k,:),ymesh(k,:),zmesh(k,:));
end
end


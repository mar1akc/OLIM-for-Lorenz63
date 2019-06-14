function handle = StrangeAttractorMesh()
global lo rho sigma beta handle colormesh

% If you want to change rho, you need to find appropriate y1 and y2 and
% replace then in the lines right below %% hit saddle

handle = [];
rho = 24.4;
sigma = 10;
beta = 8/3;
ye = [sqrt(beta*(rho-1));sqrt(beta*(rho-1));rho-1]; % equilibrium
ye1 = [-sqrt(beta*(rho-1));-sqrt(beta*(rho-1));rho-1]; % equilibrium

colormesh = [1,0,0; %  red
       0,0.8,0; % dark green
       0,0,1; %  blue;
       1,0,1; % magenta
       0.4,0.2,0; % brown
       0,0.8,0.8]; % cyan
%% hit saddle
y1 = [4.706881748620419  -0.691312411891580  27.009650570642709];
y2 = [4.706881748620418  -0.691312411891582  27.009650570642712];
n = 10;
col = cool(n);

% y0 is the set of initial conditions for trajectories passing very close
% to the saddle at the origin
g = linspace(0,1,n);
y0 = interp1([0,1],[y1;y2],g);

lo = @(t,a) [-sigma*a(1) + sigma*a(2); a(1).*(rho - a(3)) - a(2); -beta*a(3) + a(1).*a(2)];
options = odeset('AbsTol',1e-12,'RelTol',1e-12);

figure(1); clf; hold on;
grid;
view(3)
for k = 4 : 4 % the best showt is for k = 4
    [~,Y] = ode45(lo,[0,3],y0(k,:),options);
%      plot3(Y(:,1),Y(:,2),Y(:,3),'color',col(k,:),'Linewidth',2);
    fprintf('k = %d, Y(end,:) = (%d,%d,%d)\n',k,Y(end,1),Y(end,2),Y(end,3));
end

nY = sqrt(sum(Y.^2,2));
[nYmin,imin] = min(nY);

% The found trajectory passing very close to the saddle
Ysaddle = Y(1 : imin,:);

% Plot equilibria
plot3(ye(1),ye(2),ye(3),'r.','Markersize',30);
plot3(-ye(1),-ye(2),ye(3),'b.','Markersize',30);

%%
A = [-sigma,sigma,0;1,-1,-ye(1);ye(1),ye(1),-beta];
% [U,R] = LinearQpot3D(A);
U = inv(sylvester(A,A',-2*eye(3)));
[V E] = eig(U);

%% SURFACE 1 : find the border of the domain 
Nang = 72;
phi = linspace(0,2*pi,Nang + 1);
phi(end) = [];
jstar = 55;
v1 = V(:,2)*cos(phi(jstar)) + V(:,3)*sin(phi(jstar));
v1 = v1/norm(v1);
v2 = V(:,1);
a = cross(v1,v2);
%a(1:2) = -a(1:2);
events = @(t,y)MyCross(y,-a,ye);
options1 = odeset('AbsTol',1e-12,'RelTol',1e-12,'refine',10,'Events',events);

ysmall = 1e-2;
y0 = [0,-ysmall,0];
events =  @(t,y)MyCross(y,[-1,0,0]',[0,0,0]');
[~,Y0] = ode45(lo,[0,100],y0,options1);
events = @(t,y)MyCross(y,-a,[-ye(1:2);ye(3)]);
[~,Y1] = ode45(lo,[0,100],Y0(end,:),options1);
Y0(end,:) = [];
Y0 = [Y0;Y1];

% keep this for the loop
[~,Y0loop] = ode45(lo,[0,100],Y0(end,:),options1);

ha = plot3(Y0loop(:,1),Y0loop(:,2),Y0loop(:,3),'Linewidth',5,'color',colormesh(5,:));
handle = [handle,ha];
ha = plot3(-Y0loop(:,1),-Y0loop(:,2),Y0loop(:,3),'Linewidth',5,'color',colormesh(6,:));
handle = [handle,ha];

xlabel('x_1');
ylabel('x_2');
zlabel('x_3');
daspect([1,1,1])
set(gca,'FontSize',20);


ysmall = 1e-2;
y0 = [0,ysmall,0];
events = @(t,y)MyCross(y,-a,ye);
options1 = odeset('AbsTol',1e-12,'RelTol',1e-12,'refine',10,'Events',events);
[~,Y0] = ode45(lo,[0,100],y0,options1);
ha = plot3(Y0(:,1),Y0(:,2),Y0(:,3),'Linewidth',5,'color',colormesh(1,:));
handle = [handle,ha];

ha = plot3(-Y0(:,1),-Y0(:,2),Y0(:,3),'Linewidth',5,'color',colormesh(3,:));
handle = [handle,ha];

%% 
Yout = [Ysaddle;Y0];
Yin = Y0loop;
Nstart = 20;
Ncut = 72;

fname = sprintf('LorenzAttractor%d.mat',19);
data = load(fname);
b1 = data.v1;
b2 = data.v2;
r = data.x1;
z = data.y1;
ll = Y0loop(1,:) - ye';
rmin = ll*b1;
ls = Ysaddle(1,:) - ye';
rmax = ls*b1;
ind = find(r >= rmin & r <= rmax);
x = ye*ones(size(r(ind))) + b1*r(ind) + b2*z(ind);
x = [Y0loop(1,:)',x,Ysaddle(1,:)'];
rr = [rmin,r(ind),rmax]';
ri = linspace(rmin,rmax,Nstart)';

ystart = interp1(rr,x',ri);

% plot3(ystart(:,1),ystart(:,2),ystart(:,3),'.','Markersize',30);


[xmesh0,ymesh0,zmesh0] = MakeMesh(Yin,Yout,Nstart,Ncut,ystart,1);

%% SURFACE 2 : find the border of the domain 
Nang = 72;
phi = linspace(0,2*pi,Nang + 1);
phi(end) = [];
jstar = 19;
v1 = V(:,2)*cos(phi(jstar)) + V(:,3)*sin(phi(jstar));
v1 = v1/norm(v1);
v2 = V(:,1);
a = cross(v1,v2);
a(1:2) = -a(1:2);
events = @(t,y)MyCross(y,-a,[-ye(1:2);ye(3)]);
options1 = odeset('AbsTol',1e-12,'RelTol',1e-12,'refine',10,'Events',events);

y0 = Yout(end,:);
events =  @(t,y)MyCross(y,[-1,0,0]',[0,0,0]');
options1 = odeset('AbsTol',1e-12,'RelTol',1e-12,'refine',10,'Events',events);
[~,Y0] = ode45(lo,[0,100],y0,options1);
events = @(t,y)MyCross(y,-a,[-ye(1:2);ye(3)]);
options1 = odeset('AbsTol',1e-12,'RelTol',1e-12,'refine',10,'Events',events);
[~,Y1] = ode45(lo,[0,100],Y0(end,:),options1);
Y0(end,:) = [];
Y0 = [Y0;Y1];
ha = plot3(Y0(:,1),Y0(:,2),Y0(:,3),'Linewidth',5,'color',colormesh(2,:));
handle = [handle,ha];

ha = plot3(-Y0(:,1),-Y0(:,2),Y0(:,3),'Linewidth',5,'color',colormesh(4,:));
handle = [handle,ha];


% keep this for the loop
% [~,Y0loop] = ode45(lo,[0,100],Y0(end,:),options1);
% 
% plot3(Y0loop(:,1),Y0loop(:,2),Y0loop(:,3),'Linewidth',5,'color','m');
% plot3(Y0loop(1,1),Y0loop(1,2),Y0loop(1,3),'.','Markersize',40,'color','m');



ysmall = 1e-2;
y0 = [0,-ysmall,0];
events = @(t,y)MyCross(y,-a,[-ye(1:2);ye(3)]);
options1 = odeset('AbsTol',1e-12,'RelTol',1e-12,'refine',10,'Events',events);
[~,You] = ode45(lo,[0,100],y0,options1);
% plot3(You(:,1),You(:,2),You(:,3),'Linewidth',5,'color','m');

%% 
Yout = [Ysaddle;You];
Yin = Y0;
Nstart = 20;
Ncut = 72;

y = [xmesh0(:,end),ymesh0(:,end),zmesh0(:,end)];
y = y(end : -1 : 1,:);
r = (y - ones(size(y,1),1)*ye')*b1;
ind = find(r >= rmax);
yy = y(ind,:);
[yy,leny,ly] = arclength(yy);
ystart = interp1(leny,yy,linspace(0,ly,Nstart)');
ystart(1,:) = [xmesh0(end,1),ymesh0(end,1),zmesh0(end,1)];


[xmesh1,ymesh1,zmesh1] = MakeMesh(Yin,Yout,Nstart,Ncut,ystart,2);


fname = sprintf('LorenzLimitCycle_rho_%.2f.mat',rho);
data = load(fname);
gam = data.Y2;

plot3(gam(:,1),gam(:,2),gam(:,3),'Linewidth',2,'color','r');
plot3(-gam(:,1),-gam(:,2),gam(:,3),'Linewidth',2,'color','b');

end

%%
%%
function [position,isterminal,direction] = MyCross(y,a,y0)
position = (y - y0)'*a; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 1;   
end

%%

function [xmesh,ymesh,zmesh] = MakeMesh(Y0loop,Y1loop,Nstart,Ncut,ystart,snumber)
% snumber is the surface number
global lo rho sigma beta handle
[Y0loop,l0loop,l0] = arclength(Y0loop);
[Y1loop,l1loop,l1] = arclength(Y1loop);
fprintf('length(Y0loop) = %d, length(Y1loop) = %d\n',l0,l1);

gstart = linspace(0,1,Nstart)';
gcut = linspace(0,1,Ncut)';
Y0cutarray = interp1(l0loop,Y0loop,gcut*l0)';
xmesh = zeros(Nstart,Ncut);
xmesh(1,:) = Y0cutarray(1,:);
ymesh = zeros(Nstart,Ncut);
ymesh(1,:) = Y0cutarray(2,:);
zmesh = zeros(Nstart,Ncut);
zmesh(1,:) = Y0cutarray(3,:);

Y1cutarray = interp1(l1loop,Y1loop,gcut*l1)';
xmesh(end,:) = Y1cutarray(1,:);
ymesh(end,:) = Y1cutarray(2,:);
zmesh(end,:) = Y1cutarray(3,:);

xmesh(:,1) = ystart(:,1);
ymesh(:,1) = ystart(:,2);
zmesh(:,1) = ystart(:,3);
for k = 2 : Ncut    
%     fprintf('Cut # %i\n',k);
    Y0cut = Y0cutarray(:,k);
    Y1cut = Y1cutarray(:,k);
    bcut = lo(0,Y0cut);
    bcut = bcut/norm(bcut);
    a = Y1cut - Y0cut;
    a = a/norm(a);
    aa = cross(a,bcut);
    aa = aa/norm(aa);
    b = cross(aa,a);
    events = @(t,y)MyCross(y,b,Y0cut);
    options1 = odeset('AbsTol',1e-12,'RelTol',1e-12,'refine',10,'Events',events);
    yend = zeros(Nstart,3);
    yend(1,:) = Y0cut';
    yend(end,:) = Y1cut';
    for j = 2 : Nstart - 1
        y0 = ystart(j,:);
        [~,Y] = ode45(lo,[0,100],y0,options1);
        aux = Y(end,:);
        yend(j,:) = Y(end,:);
    end
    [yend,lyend,ly] = arclength(yend);
    ystart = interp1(lyend,yend,linspace(0,ly,Nstart)');
    xmesh(:,k) = ystart(:,1);
    ymesh(:,k) = ystart(:,2);
    zmesh(:,k) = ystart(:,3);
end
plotmesh(xmesh,ymesh,zmesh,snumber);
plotmesh(-xmesh,-ymesh,zmesh,snumber + 2);
end
%%
function [Y0,lxyz,lY] = arclength(Y0)
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
function plotmesh(xmesh,ymesh,zmesh,cind)
global handle colormesh
[Nstart,Ncut] = size(xmesh);
for k = 1 : Ncut
    ha = plot3(xmesh(:,k),ymesh(:,k),zmesh(:,k),'color',colormesh(cind,:));
    handle = [handle,ha];
end
for k = 1 : Nstart
    ha = plot3(xmesh(k,:),ymesh(k,:),zmesh(k,:),'color',colormesh(cind,:));
    handle = [handle,ha];
end
end

% %%
% function [u,v,w] = RotateSpace(xmesh,ymesh,zmesh)
% [Nstart,Ncut] = size(xmesh);
% xyz = [xmesh(:),ymesh(:),zmesh(:)];
% co = cos(pi/4); si = sin(pi/4);
% M = [co,si,0;-si,co,0;0,0,1];
% uvw = (M*xyz')';
% u = reshape(uvw(:,1),Nstart,Ncut);
% v = reshape(uvw(:,2),Nstart,Ncut);
% w = reshape(uvw(:,3),Nstart,Ncut);
% end



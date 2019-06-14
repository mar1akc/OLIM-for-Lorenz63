function MyPolarPlot() %xmesh,ymesh,zmesh,Nrad,Nloop,u)
rho = 15.0;
sigma = 10;
beta = 8/3;

ye = [sqrt(beta*(rho-1));sqrt(beta*(rho-1));rho-1]; % equilibrium

fname = sprintf('Qpot2Drho%.2f.txt',rho);
u = load(fname);
[Nrad,Nloop] = size(u);
u = [u,u(:,1)];
Nloop = Nloop + 1;

fname = sprintf('Mesh2Drho%.2f.txt',rho);
mesh = load(fname);
xmesh = mesh(1 : Nrad,:);
ymesh = mesh(Nrad + 1 : 2*Nrad,:);
zmesh = mesh(2*Nrad + 1 : end,:);
xmesh = [xmesh,xmesh(:,1)];
ymesh = [ymesh,ymesh(:,1)];
zmesh = [zmesh,zmesh(:,1)];

I = 1 : 4 : Nloop;
J = 1 : 4 : Nrad;

u1 = u(J,I);
Nloop = length(I);
Nrad = length(J);

ind = find( u > 1e5);
u(ind) = NaN;
t = linspace(0,2*pi,Nloop)'; % angle

Rmax = sqrt(sum(([xmesh(end,I);ymesh(end,I);zmesh(end,I)] - ye*ones(1,Nloop)).^2,1))';
for i = 1 : Nloop
    r = linspace(0,Rmax(i),Nrad)';
    x(:,i) = r*cos(t(i));
    y(:,i) = r*sin(t(i));
    qcyc(i) = u(end,i);
end

fprintf('U at the limit cycle: mean = %d, mi = %d, max = %d\n',mean(qcyc),min(qcyc),max(qcyc));

figure(3); clf; hold on; grid;
surf(x',y',u1','Edgecolor','none');
view(2)
plot(0,0,'.','color','r','Markersize',30);
plot(Rmax.*cos(t),Rmax.*sin(t),'r','Linewidth',3);
daspect([1,1,1])
xlabel('v1','Fontsize',20);
ylabel('v2','Fontsize',20);
colorbar;
set(gca,'Fontsize',20);
end

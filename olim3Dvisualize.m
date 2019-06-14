function handle = olim3Dvisualize()
close all
load_u = 1;
rho = 24.4;
%load file with parameters
fname = sprintf('parameters_rho%.2f.txt',rho);
param = load(fname);
XMIN = param(1);
XMAX = param(2);
YMIN = param(3);
YMAX = param(4);
ZMIN = param(5);
ZMAX = param(6);
NX = param(7);
NY = param(8);
NZ = param(9);
sigma = param(10);
beta = param(11);
rho = param(12);
ye = [param(13),param(14),param(15)];
K = param(16);
subsample = param(17);
rmax = norm(ye);

handle = [];


NX = (NX - 1)/subsample + 1;
NY = (NY - 1)/subsample + 1;
NZ = (NZ - 1)/subsample + 1;

if load_u == 1 
    fname = sprintf('LorenzQpot_rho%.2f.txt',rho);
    % load the file with the quasi-potential
    uvec = load(fname);
    ind = find(uvec > 1e5);
    uvec(ind) = NaN;
    umax = max(uvec);
    u1 = reshape(uvec,NY*NX,NZ);
    for k = 1 : NZ
        u(:,:,k) = reshape(u1(:,k),[NY,NX]);
    end
end
figure(1); hold on; grid;
xx = linspace(XMIN,XMAX,NX);
yy = linspace(YMIN,YMAX,NY);
zz = linspace(ZMIN,ZMAX,NZ);
hz = zz(2)-zz(1);
[x,y,z] = meshgrid(xx,yy,zz);

ha = plot3(ye(1),ye(2),ye(3),'r.','Markersize',30);
handle = [handle,ha];
ha = plot3(-ye(1),-ye(2),ye(3),'r.','Markersize',30);
handle = [handle,ha];

% isoval = linspace(0,0.97*umax,nval + 1);
% isoval(1) = [];
if rho < 1
    isoval = [20, 40];
    nval = 2;
end
if rho > 1 & rho < 13.926
    uorigin = interp3(x,y,z,u,0,0,0);
    fprintf('rho = %d, uorigin = %d\n',rho,uorigin);
    isoval(1) = 10;
    isoval(2) = uorigin-0.05;
    nval = 2;
end
col = ['b','r'];

if (rho > 13.926 & rho < 24.74 ) 
    [loop,llop,len] = find_saddle_cycle(rho);
    qloop = interp3(x,y,z,u,loop(:,1),loop(:,2),loop(:,3));
    maxqloop = max(qloop);
    minqloop = min(qloop);
    meanqloop = mean(qloop);
    fprintf('limit cycle: qmax = %d, qmin = %d, qmean = %d\n',maxqloop,minqloop,meanqloop);
    if abs(rho - 15) < 1e-14        
        isoval = [8,minqloop - 0.01, 20];
        col = ['g','b','r'];
    end
    if abs(rho - 20) < 1e-14
        isoval = [minqloop/2,minqloop - 0.01];
        col = ['b','r'];
    end

    if abs(rho - 24.4) < 1e-14
        isoval = [minqloop - 0.01,2,20];
         col = ['g','b','r'];
    end
    nval = length(isoval);
    ha = plot3(loop(:,1),loop(:,2),loop(:,3),'r','Linewidth',3);
    handle = [handle,ha];
    ha = plot3(-loop(:,1),-loop(:,2),loop(:,3),'r','Linewidth',3);
    handle = [handle,ha];

end
    
if abs(rho - 100.75) < 1e-15
    fname = sprintf('LorenzLimitCycle_rho_%.2f.mat',rho);
    data = load(fname);
    loop = data.Y2;
    qloop = interp3(x,y,z,u,loop(:,1),loop(:,2),loop(:,3));
    maxqloop = max(qloop);
    minqloop = min(qloop);
    meanqloop = mean(qloop);
    fprintf('limit cycle: qmax = %d, qmin = %d, qmean = %d\n',maxqloop,minqloop,meanqloop);
    nval = 2;
    isoval = [4,8];
    col = ['b','r'];
    ha = plot3(loop(:,1),loop(:,2),loop(:,3),'r','Linewidth',3);
    handle = [handle,ha];

end
    
kstep = 20;
for i = 1 : nval
    figure(1)
    p = patch(isosurface(x,y,z,u,isoval(i)));
    p.FaceColor = col(i);
    p.EdgeColor = 'none';
    alpha(0.2)
    fprintf('isoval(%d) = %d\n',i,isoval(i))
    [faces,verts] = isosurface(x,y,z,u,isoval(i));
    v3 = verts(:,3);
    vmax = max(v3);
    vmin = min(v3);
    kmin = ceil((vmin-ZMIN)/hz);
    kmax = floor((vmax -ZMIN)/hz);
    handle = [handle,p];
    if rho < 90
        for k = kmin : kstep : kmax
            figure(2); hold on;
            [c,hc] = contour(xx,yy,u(:,:,k),[isoval(i),isoval(i)]);
            if ~isempty(c)
                m = size(c,2);
                icur = 1;
                while icur < m
                    figure(1);
                    idat = c(2,icur);
                    ha = plot3(c(1,icur + 1 : icur + idat),c(2,icur + 1 : icur + idat),zz(k)*ones(1,idat),...
                        'color',col(i),'Linewidth',1);
                    icur = icur + idat + 1;
                    handle = [handle,ha];
                end
                drawnow;
            end
        end
    end
        
end



daspect([1,1,1])
set(gca,'FontSize',20);
xlabel('x_1');
ylabel('x_2');
zlabel('x_3');
view(3);
if rho > 13 & rho < 24.06
    [zmax,imax] = max(loop(:,3));
    v = loop(imax,3) - ye;
    v = v/norm(v);
    lo = @(t,a) [-sigma*a(1) + sigma*a(2); a(1).*(rho - a(3)) - a(2); -beta*a(3) + a(1).*a(2)];
    myevents = @(t,y)ReachEq(y,ye);
    options = odeset('AbsTol',1e-12,'RelTol',1e-12,'Events',myevents);
    yi = loop(imax,:) - 1e-3*v;
    fprintf('initial point: %d\t%d\t%d\n',yi(1),yi(2),yi(3));
    [~,Y] = ode45(lo,[0,100],yi,options);
    ha = plot3([Y(:,1)],[Y(:,2)],[Y(:,3)],'color',[0,0,0.6],'Linewidth',2);
    handle = [handle,ha];
    yo = loop(imax,:) + 1e-3*v;
    [~,Y] = ode45(lo,[0,100],yo,options);
    ha = plot3([Y(:,1)],[Y(:,2)],[Y(:,3)],'color',[0,0,0.6],'Linewidth',2);
    handle = [handle,ha];
    fprintf('yi: %d\t%d\t%d\n',yi(1),yi(2),yi(3))
    fprintf('yo: %d\t%d\t%d\n',yo(1),yo(2),yo(3))
    
    % plot a MAP arriving at saddle
    fname = sprintf('MAP_rho%.2f.txt',rho);
    map = load(fname);
    ha = plot3(map(:,1),map(:,2),map(:,3),'color',[0.6,0,0],'Linewidth',2);
    handle = [handle,ha];

end    
if rho > 1 & rho < 13
    % plot trajectory emanating from saddle
    lo = @(t,a) [-sigma*a(1) + sigma*a(2); a(1).*(rho - a(3)) - a(2); -beta*a(3) + a(1).*a(2)];
    myevents = @(t,y)ReachEq(y,ye);
    options = odeset('AbsTol',1e-12,'RelTol',1e-12,'Events',myevents);
    [~,Y] = ode45(lo,[0,100],[0,0.01,0],options);
    ha = plot3([0;Y(:,1)],[0;Y(:,2)],[0;Y(:,3)],'color',[0,0,0.6],'Linewidth',3);
    handle = [handle,ha];

    % plot a MAP arriving at saddle
    fname = sprintf('MAP_rho%.2f.txt',rho);
    map = load(fname);
    ha = plot3(map(:,1),map(:,2),map(:,3),'color',[0.6,0,0],'Linewidth',3);
    handle = [handle,ha];
end
if rho < 1
    figure(1);
    hold on;
    Nang = 72;
    t = linspace(0,2*pi,Nang + 1);
    t(end) = [];
    r = 0.01 : 0.02 : rmax;
    nr = length(r);
    s2 = 1/sqrt(2);
    rad = 0.1;
    rad2 = rad*rad;
    for k = 1 : Nang
        for i = nr : -1 : 1
            c1 = r(i)*cos(t(k));
            c2 = r(i)*sin(t(k));
            ind = find((((verts(:,1) + verts(:,2))*s2 - c1).^2 + (verts(:,3) - c2).^2) < rad2);
            if length(ind) > 20
                R(k) = r(i);
                Zeta(k) = mean((verts(ind,1) - verts(ind,2))*s2);
                break;
            end
        end
    end

    % plot trajectories and MAPs
    lo = @(t,a) [-sigma*a(1) + sigma*a(2); a(1).*(rho - a(3)) - a(2); -beta*a(3) + a(1).*a(2)];
    options = odeset('AbsTol',1e-12,'RelTol',1e-12,'Events',@events);
    a1 = s2*(R.*cos(t) + Zeta) + ye(1);
    a2 = s2*(R.*cos(t) - Zeta) + ye(2);
    a3 = R.*sin(t) + ye(3);
    ha = plot3(a1,a2,a3,'.','Markersize',20);
    handle = [handle,ha];
    for k = 1 :  Nang
        [~,Y] = ode45(lo,[0,100],[a1(k),a2(k),a3(k)],options);
        plot3(s2*(Y(:,1)+Y(:,2)),Y(:,3),s2*(Y(:,1)-Y(:,2)));
        ha = plot3(Y(:,1),Y(:,2),Y(:,3),'color',[0,0,0.4],'Linewidth',1);
        handle = [handle,ha];
        MAP = gmam_lorenz(ye',[a1(k);a2(k);a3(k)],sigma,beta,rho);
        ha = plot3(MAP(:,1),MAP(:,2),MAP(:,3),'color',[0.4,0,0],'Linewidth',1);
        handle = [handle,ha];
    end
end
end
%%
%%
function [position,isterminal,direction] = events(t,y)
position = norm(y) - 1e-3; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 0;   
end

%%
function [position,isterminal,direction] = ReachEq(y,ye)
position = min([norm(y - ye') - 1e-3,norm(y - [-ye(1);-ye(2);ye(3)]) - 1e-3]); % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 0;   
end

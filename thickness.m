function w = thickness(y0)
rho = 24.4;
sigma = 10;
beta = 8/3;
ye = [sqrt(beta*(rho-1));sqrt(beta*(rho-1));rho-1]; % equilibrium
ye1 = [-sqrt(beta*(rho-1));-sqrt(beta*(rho-1));rho-1]; % equilibrium

lo = @(t,a) [-sigma*a(1) + sigma*a(2); a(1)*(rho - a(3)) - a(2); -beta*a(3) + a(1)*a(2)];

% Some points lying on the strange attractor
% y0 = [ -10.165034795972931 -15.462661714122202  17.969596907054004]';
% y0 = [ 11.010571987772378  13.555096197136503  23.937827067336414]';
% y0 = [13.086843421346815  10.040427061932833  32.948498131434299]';
% y0 = [ 0.396509769359544   0.706269596790777  10.842731108830900]';
% y0 = [ -6.148951223738913  -7.447219547641666  19.156567468262569]';

% if flag == 1 draw a figure, otherwise do not
flag = 1;

b0 = lo(0,y0);
b0 = b0/norm(b0);
v1 = [1;0;0];
v1 = v1 - (v1'*b0)*b0;
v1 = v1/norm(v1);
v2 = cross(b0,v1);


events = @(t,y)mycross(y,b0,y0);
options = odeset('AbsTol',1e-8,'RelTol',1e-8,'Events',events);

flag = 1;

if flag == 1
    fprintf('Computing the long trajectory ... ');
    [T,Y,Ti,Yi,I] = ode45(lo,[0,10000],[0,1e-2,0],options);
    fprintf('done\n');
    save('thickdata.mat','Yi','y0','b0','v1','v2');
else
    data = load('thickdata.mat');
    Yi = data.Yi;
    y0 = data.y0;
    v1 = data.v1;
    v2 = data.v2;
end
c1 = (Yi - ones(size(Yi,1),1)*y0')*v1;
c2 = (Yi - ones(size(Yi,1),1)*y0')*v2;
r = 0.25;
xp = -r : 0.01*r : r;
ind = find(abs(c1) < r & abs(c2) < r);
d1 = c1(ind);
d2 = c2(ind);
p = polyfit(d1,d2,1);
i1 = find(d2 > polyval(p,d1));
p1 = polyfit(d1(i1),d2(i1),1);
i2 = find(d2 < polyval(p,d1));
p2 = polyfit(d1(i2),d2(i2),1);
pp = [-1/p1(1),0];
if abs(p1(2)) < abs(p2(2))
    xc = -p2(2)/(p1(1) + 1/p1(1));
else
    xc = -p1(2)/(p1(1) + 1/p1(1));
end
yc = polyval(pp,xc);
w = norm([xc,yc]);
fprintf('thickness = %d\n',w);

if flag == 1
    figure;
    hold on;
    grid;
    plot(c1,c2,'.');
    axis([-r,r,-r,r]);
    plot(xp,polyval(p,xp))
    plot(d1(i1),d2(i1),'r.')
    plot(xp,polyval(p1,xp))
    plot(xp,polyval(p2,xp));
    plot(xp,polyval(pp,xp));
    daspect([1,1,1]);
end
end
                
%%
function [position,isterminal,direction] = mycross(y,a,y0)
position = (y - y0)'*a; % The value that we want to be zero
isterminal = 0;  % Halt integration 
direction = 1;   
end





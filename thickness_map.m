function thickness_map()
% StrangeAttractorMesh();
rho = 24.4;
sigma = 10;
beta = 8/3;
ye = [sqrt(beta*(rho-1));sqrt(beta*(rho-1));rho-1]; % equilibrium
ye1 = [-sqrt(beta*(rho-1));-sqrt(beta*(rho-1));rho-1]; % equilibrium
lo = @(t,a) [-sigma*a(1) + sigma*a(2); a(1)*(rho - a(3)) - a(2); -beta*a(3) + a(1)*a(2)];
options = odeset('AbsTol',1e-8,'RelTol',1e-8);
Tmax = 1000;
y0 = [0,1e-2,0]; % [1.912041e+00    1.411822e+00    1.803113e+01];
[T,Y] = ode45(lo,[0,Tmax],y0,options);
t = 10;
k = 1;
col = jet(10);
while t <= Tmax
    [dmin,imin] = min(abs(T-t));
    y0 = Y(imin,:)';
    w(k) = thickness(y0);
    lw = max(1,min(floor(-log(w(k))/log(10)),10));
    figure(1);
    plot3(y0(1),y0(2),y0(3),'.','Markersize',30,'color',col(lw,:));
    drawnow;
    record_thickness(y0,w(k),lw);
    t = t + 5;
end
end
%%
function record_thickness(y0,w,lw)
fid = fopen('thickness24p4.txt','a');
fprintf(fid,'%d\t%d\t%d\t%d\t%d\n',y0(1),y0(2),y0(3),w,lw);
fclose(fid);
end

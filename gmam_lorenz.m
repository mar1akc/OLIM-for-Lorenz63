function MAP = gmam_lorenz(xi,xf,sigma,beta,rho)
% Heymann's and V.E.'s geometric minimum action method
% compites the minimum action path from
n = 100;


    n1 = n-1;
    n2 = n - 2;
    % define the initial
    t = linspace(0,1,n);
    x = interp1([0 1], [xi xf]',t)';
    [x,l] = reparametrization(x); % Note: |x'| = l;

    % Set up figure
    % figure(1);
    % clf; hold on; grid;
%     mypath = plot3(x(1,:),x(2,:),x(3,:),'k');
%     plot3([xi(1),xs(1),xf(1)],[xi(2),xs(2),xf(2)],[xi(3),xs(3),xf(3)],'.','Markersize',30);
    view(3);

    % Parameters of the method
    nor=10;
    % parameters in gMAM
    tau = 5.0e-3; % time step
    tol = 1.0e-6;
    %% start
    k=0;
    h = (1/n1); % the path is parametrized from 0 to 1. h = step in parameter
    r = tau*n1*n1*ones(n,1);
    D2 = spdiags([r -2*r r],-1:1,n2,n2);
    while nor > tol
        xold = x;
        dxa = 0.5*(circshift(x,[0,-1])-circshift(x,[0,1]))/h;   % x' along the path
        [b1,b2,b3] = bfield(x,sigma,beta,rho);
        lam = sqrt(b1.^2 + b2.^2 + b3.^2)/l; % |b|/|x'|
        dlam = 0.5*(circshift(lam,[0,-1])-circshift(lam,[0,1]))/h;
        [b1x,b1y,b1z,b2x,b2y,b2z,b3x,b3y,b3z] = bgrad(x,sigma,beta,rho);
        dbx = (b2x - b1y).*dxa(2,:) + (b3x - b1z).*dxa(3,:); % [[(grad b)^T - grad b ]*x'] -- x-component
        dby = (b1y - b2x).*dxa(1,:) + (b3y - b2z).*dxa(3,:); % [[(grad b)^T - grad b ]*x'] - y-component
        dbz = (b1z - b3x).*dxa(1,:) + (b2z - b3y).*dxa(2,:); % [[(grad b)^T - grad b ]*x'] - z-component
        bx = b1x.*b1 + b2x.*b2 + b3x.*b3; % (grad b)^T b -- x-component
        by = b1y.*b1 + b2y.*b2 + b3y.*b3; % (grad b)^T b -- y-component
        bz = b1z.*b1 + b2z.*b2 + b3z.*b3; % (grad b)^T b -- z-component
        % linearly implicit scheme
        mymatr=(eye(n2) - diag(lam(2:n1).^2)*D2);

        rhsx = x(1,2:n1) + tau*(lam(2:n1).*dbx(2:n1) - bx(2:n1) + lam(2:n1).*dlam(2:n1).*dxa(1,2:n1));
        rhsy = x(2,2:n1) + tau*(lam(2:n1).*dby(2:n1) - by(2:n1) + lam(2:n1).*dlam(2:n1).*dxa(2,2:n1));
        rhsz = x(3,2:n1) + tau*(lam(2:n1).*dbz(2:n1) - bz(2:n1) + lam(2:n1).*dlam(2:n1).*dxa(3,2:n1));

        rhsx(1) = rhsx(1)+tau*x(1,1)*(n1*lam(2))^2;
        rhsy(1) = rhsy(1)+tau*x(2,1)*(n1*lam(2))^2;
        rhsz(1) = rhsz(1)+tau*x(3,1)*(n1*lam(2))^2;
        rhsx(n2) = rhsx(n2)+tau*x(1,end)*(n1*lam(n1))^2;
        rhsy(n2) = rhsy(n2)+tau*x(2,end)*(n1*lam(n1))^2;
        rhsz(n2) = rhsz(n2)+tau*x(3,end)*(n1*lam(n1))^2;

        x(1,2:n1) = (mymatr\rhsx')';
        x(2,2:n1) = (mymatr\rhsy')';
        x(3,2:n1) = (mymatr\rhsz')';

        [x,l] = reparametrization(x);
        k = k + 1;
        nor = norm(x - xold)/tau;
%         fprintf('iter # %d:  res = %d\n',k,nor);

    %     figure(1);
%         set(mypath,'Xdata',x(1,:),'Ydata',x(2,:),'Zdata',x(3,:));
%         drawnow;
    end
%     set(mypath,'Xdata',x(1,:),'Ydata',x(2,:),'Zdata',x(3,:));
%     drawnow;
    MAP = x';
%     fname = sprintf('MAP%d.mat',ipath);
%     save(fname,'x');

%     [~,isort] = sort(sqrt(sum(x.^2,1)),'ascend');
%     [~,Y] = ode45(@GS,[0,1000],x(:,isort(1))');
%     plot3(Y(:,1),Y(:,2),Y(:,3),'b','Linewidth',2);
%     [~,Y] = ode45(@GS,[0,1000],x(:,isort(2))');
%     plot3(Y(:,1),Y(:,2),Y(:,3),'b','Linewidth',2);
    % compute quasi-potential
    dxa = (circshift(x,[0,-1])-x)/h;   % x' along the path
    [b1,b2,b3] = bfield(x,sigma,beta,rho);
    qp=0;
    for i = 1 : n1
        bb= 0.5*([b1(i),b2(i),b3(i)] + [b1(i + 1),b2(i + 1),b3(i + 1)]);
        qp=qp+(norm(bb)*norm(dxa(:,i))-bb*dxa(:,i))*h;
    end
    fprintf('quasi-potential at xf is %.4e\n',qp); 
end
%%
function [x,l] = reparametrization(x)
% x in an d by n array. Row i of x is the coordinate x_i along the path
% returns a uniformly reparametrized path and its length l
t = linspace(0,1,size(x,2));
dx = x - circshift(x,[0,1]);
dx(:,1) = zeros(size(x,1),1);
lxy = cumsum(sqrt(sum(dx.^2,1)));
l = lxy(end);
x = interp1(lxy/l,x',t)';
end
%%
function [b1 b2 b3] = bfield(x,sigma,beta,rho)
% Input is an d by n array 
b1 = sigma*(x(2,:) - x(1,:));
b2 = x(1,:).*(rho - x(3,:)) - x(2,:);
b3 = x(1,:).*x(2,:) - beta*x(3,:);
end
%%
function [b1x b1y b1z b2x b2y b2z b3x b3y b3z] = bgrad(x,sigma,beta,rho)
% Input is an d by n array 
b1x = -sigma*ones(1,size(x,2)); b1y = sigma*ones(1,size(x,2)); b1z = zeros(1,size(x,2));
b2x = rho*ones(1,size(x,2)) - x(3,:); b2y = -ones(1,size(x,2)); b2z = -x(1,:);
b3x = x(2,:); b3y = x(1,:); b3z = -beta*ones(1,size(x,2));
end


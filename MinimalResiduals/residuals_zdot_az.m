%% Interpolant and residual curves
% Yalcin Kaya 11 Jan 2014, graphs formatted further on 19-20 Nov 2017
% zdot = a*z

clear all

fsize = 60; msize = 40;
a = 3;
numdiv = 50;
scale_flag = 0;

f = @(t,y) a*y(1);
tol = 1e-8;
opts = odeset('RelTol',tol,'AbsTol',tol);
sol = ode15s(f, [0,1], 1, opts);
%sol = ode45(f, [0,1], 1, opts);
%sol = ode113(f, [0,1], 1, opts);
tr = RefineMesh( sol.x, numdiv );
n = length(tr);
res = zeros(1,n);
[y,dy] = deval(sol,tr);
for i=1:n
    res(:,i) = dy(:,i) - f(tr(i),y(:,i));
end
figure(1); clf; plot( tr, y, 'k.','MarkerSize', msize ), grid on, 
xlabel('$t$','Interpreter','latex','Fontsize',40), 
ylabel('{$x(t)$}','Interpreter','latex','Fontsize',40);
set(gca,'FontSize', fsize);
figure(2); clf; plot( tr, res(1,:), 'b.','MarkerSize', msize ), grid on
title(sol.solver,'Interpreter','latex','Fontsize',80);
xlabel('$t$','Interpreter','latex','Fontsize',80), 
ylabel('$u_M(t)$','Interpreter','latex','Fontsize',80);
set(gca,'FontSize', fsize);

tnode = sol.x;
znode = sol.y;

%% L2-minimization of the residual

t = zeros(length(tnode)-1,numdiv+1);
x = zeros(length(tnode)-1,numdiv+1);
u = zeros(length(tnode)-1,numdiv+1);

figure(3), clf, hold on, grid on

for i = 1:length(tnode)-1
    h = tnode(i+1)-tnode(i);
    deltat = h/numdiv;
    for j = 1:numdiv+1
        t(i,j) = tnode(i) + (j-1)*deltat;
        u(i,j) = -2*a*(znode(i+1) - exp(a*h)*znode(i)) ...
            /(1 - exp(2*a*h)) * exp(a*(tnode(i+1)-t(i,j)));
        x(i,j) = -(znode(i+1) - exp(-a*h)*znode(i)) ...
            /(exp(-a*h) - exp(a*h)) * exp(a*(t(i,j) - tnode(i))) ...
            + u(i,j)/(-2*a);
    end
    figure(3), plot( t(i,:), u(i,:), 'b.','MarkerSize', msize ), 
    title(sol.solver,'Interpreter','latex','Fontsize',80);
    xlabel('$t$','Interpreter','latex','Fontsize',80), 
    ylabel('$u_{L^2}(t)$','Interpreter','latex','Fontsize',80);
    set(gca,'FontSize', fsize); 
    u_L2norm(i) = sqrt(-2*a*(znode(i+1) - exp(a*h)*znode(i))^2 ...
        /(1 - exp(2*a*h)));
end

%% Linfinity-minimization of the residual

figure(4), clf, hold on, grid on,

for i = 1:length(tnode)-1
    h = tnode(i+1)-tnode(i);
    deltat = h/numdiv;
    for j = 1:numdiv+1
%         t(i,j) = tnode(i) + (j-1)*deltat;
        u(i,j) = -a*(znode(i+1) - exp(a*h)*znode(i)) / (1 - exp(a*h));
        x(i,j) = (znode(i) + u(i,j)/a) * exp(a*(t(i,j) - tnode(i))) ...
            - u(i,j)/a;
    end
    figure(4), plot( t(i,:), u(i,:), 'b.','MarkerSize', msize ),
    title(sol.solver,'Interpreter','latex','Fontsize',50); %axis([0 1 0 2.5e-5])
    xlabel('$t$','Interpreter','latex','Fontsize',55), 
    ylabel('$u_{L^\infty}(t)$','Interpreter','latex','Fontsize',55);
    set(gca,'FontSize', fsize); 
end

if scale_flag==1
    for i=1:3
        if i==1; figure(2); end
        if i==2; figure(3); end
        if i==3; figure(4); end
        % get ylim
        yl=ylim;
        % get order of magnitude
        e=log10(yl(2));
        e=sign(e)*floor(abs(e)) - 1;
        if i==1; e=-2; end
        % get and rescale yticks
        yt=get(gca,'ytick')/10^e;
        % create tick labels
        ytl=cell(size(yt));
        for j=1:length(yt)
            % the space after the percent gives the same size to positive and
            % negative numbers. The number of decimal digits can be changed.
            if i==2 || i==3; ytl{j}=sprintf('% 1.0f',yt(j)); end
            if i==1; ytl{j}=sprintf('% 1.0f',yt(j)); end
        end
        % set tick labels
        set(gca,'yticklabel',ytl);
        % place order of magnitude
        fs = get(gca,'fontsize');
        set(gca,'units','normalized');
        xl = xlim;
        text(xl(1),yl(2),sprintf('\\times10^{%d}',e),...
            'fontsize',fs,'VerticalAlignment','bottom');
    end
end
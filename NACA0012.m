clc;
clear;
close all;

set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',20)
set(0,'DefaultLineLineWidth',1.5);
%% Airfoil Characteristics

c = 1; % Chord length
tc = 12./100; % Thickness to Chord Ratio

%% Joukowski circle 

a = (c./4) * ( 1 + 2.*tc.*4/(3.*sqrt(3)) )  .* (1 + tc.*4/(3.*sqrt(3)) )^(-2);

b = a .* tc .* 4/(3.*sqrt(3));

R = a + b;
Zeta_origin = -b; % forward facing
Zeta_origin = b; % backward facing

Zeta_circ = R*exp(1i*linspace(0,2*pi,200)) + Zeta_origin;


figure,
subplot(1,2,1), hold all
plot(complex(0),'k.') % this is (0,0)
plot(complex(Zeta_origin),'bx') % this is the center of the circle 
plot(complex(Zeta_circ),'b')
axis equal tight, axis([-1.5*R,1.5*R,-1.5*R,1.5*R]) 
grid on, grid minor
xlabel('\eta'), ylabel('\xi')


subplot(1,2,2), hold all
plot(joukowski_mapping(Zeta_circ,a),'b')
view(2), axis equal tight, axis([-0.75*c,0.75*c,-0.75*c,0.75*c]) 
grid on, grid minor
xlabel('x'), ylabel('y')

%% Flow Field

alpha = 7; % AoA in degrees
alpha = deg2rad(alpha);

Re = 2e4;
uinf = Re .* 1.56e-5 / c;

Gamma = -4 .* pi .* R .* uinf .* sin(alpha); % alpha

n = 500;
bound = 3*R;


[eta, xi] = meshgrid(linspace(-bound,bound,n), linspace(bound,-bound,n));
zeta = eta + 1i .* xi;

cond = abs(zeta-Zeta_origin)<=R;
zeta(cond) = (R).*exp(-1i.* linspace(0,2*pi,length(zeta(cond))) ) + Zeta_origin;

eta = real(zeta);
xi = imag(zeta);

x = real(joukowski_mapping(zeta,a));
y = imag(joukowski_mapping(zeta,a));

w = flow_field_gamma(alpha, uinf, R, Gamma, zeta, Zeta_origin);

psi = imag(w);

U_zeta = velocity_field_gamma(alpha, uinf, R, Gamma, zeta, Zeta_origin);

cond2 = abs((x+1i.*y) )<=1/(1.5*100);
U_zeta(cond) = 0;
U_mod(cond2) = 0;
U_real = U_zeta ./(1 - a.^2./(zeta.^2) );
Umod = abs(U_real);

% cond = abs(Umod)> mean(Umod,"all") + 2*std(Umod(:));
% Umod(cond) = 0; %mean(Umod,"all") + 2*std(Umod(:));


%%
nlevels = 6;

levels = imag(flow_field_gamma(alpha, uinf, R, Gamma, (-10*R + 1i.*linspace(0,1.5*R,nlevels)).*exp(1i.*alpha) , Zeta_origin));

if ~ismember(0,levels)
    levels = cat(2,0,levels);
end

levels = sort(cat(2,-levels(2:end),levels));

figure
hold on
colormap hot;
h1 = pcolor(eta,xi,abs(U_zeta));
contour(eta,xi,psi,levels,'k-')
plot(real(Zeta_circ),imag(Zeta_circ),'b','LineStyle','-')
hold off
axis square
grid on, grid minor
set(h1, 'EdgeColor', 'none');
colorbar
% clim([0 max(abs(U_zeta),[],"all")])
savefig('FlowContourJouk.fig')

%%

figure
hold on
colormap hot;
%h2 = pcolor(x,y,abs(U_real));
h2 = pcolor(x,y,Umod);
contour(x,y,psi,levels,'k-')
plot(real(joukowski_mapping(Zeta_circ,a)),imag(joukowski_mapping(Zeta_circ,a)),'b','LineStyle','-')
hold off
axis square
grid on, grid minor
set(h2, 'EdgeColor', 'none');
colorbar


%%


% narrows = nlevels*4;
% 
% [eta2, xi2] = meshgrid(linspace(-bound,bound,narrows), linspace(bound,-bound,narrows));
% 
% zeta2 = eta2 + 1i.*xi2;
% 
% 
% cond = abs(zeta2-Zeta_origin)<=R;
% length(cond)
% zeta2(cond) = R.*exp(-1i.* linspace(0,2*pi,length(zeta2(cond))) ) + Zeta_origin;
% 
% eta2 = real(zeta2);
% xi2 = imag(zeta2);
% 
% U_zeta = velocity_field_gamma(alpha, uinf, R, Gamma, zeta2, Zeta_origin);
% u_zeta = real(U_zeta);
% v_zeta = -imag(U_zeta);
% 
% U_real = U_zeta ./(1 - a.^2./(zeta2.^2) );
% u_real = real(U_real);
% v_real = -imag(U_real);
% 
% 
% x2 = real(joukowski_mapping(zeta2,a));
% y2 = imag(joukowski_mapping(zeta2,a));
% min(x2,[],"all")
% min(y2,[],"all")
% 
% zflow = linspace(min(x2,[],"all"),0,5) .* exp(1i.*alpha);
% 
% 
% figure
% hold on
% fill(real(Zeta_circ),imag(Zeta_circ),'b','LineStyle','none')
% quiver(eta2,xi2,u_zeta,v_zeta,'k')
% hold off
% axis square
% grid on, grid minor
% savefig('FlowArrowsJouk.fig')
% 
% cond = abs(zeta-Zeta_origin)<=R;
% zeta(cond) = R.*exp(-1i.* linspace(0,2*pi,length(zeta(cond))) ) + Zeta_origin;
% figure
% colormap hot;
% h = pcolor(real(joukowski_mapping(zeta,a)),imag(joukowski_mapping(zeta,a)),abs( velocity_field_gamma(alpha, uinf, R, Gamma, zeta, Zeta_origin) ./(1 - a.^2./(zeta.^2) ) ));
% set(h, 'EdgeColor', 'none');
% colorbar
% 
% figure
% hold on
% plot(complex(zflow),"k--")
% fill(real(joukowski_mapping(Zeta_circ,a)),imag(joukowski_mapping(Zeta_circ,a)),'b','LineStyle','none')
% quiver(x2,y2,u_real,v_real,'k')
% hold off
% legend({"freestream flow", "", ""},'Location','best')
% axis square
% grid on, grid minor



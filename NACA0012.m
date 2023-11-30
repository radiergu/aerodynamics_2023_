clc;
clear;
close all;

set(0,'DefaultLineLineWidth',1.25);
%% Airfoil Characteristics

c = 1; % Chord length
tc = 12./100; % Thickness to Chord Ratio

%% Joukowski circle 

a = (c./4) * ( 1 + 2.*tc.*4/(3.*sqrt(3)) )  .* (1 + tc.*4/(3.*sqrt(3)) )^(-2)

b = a .* tc .* 4/(3.*sqrt(3))

R = a + b;
Zeta_origin = b;
Zeta_circ = R*exp(1i*linspace(0,2*pi,200)) + Zeta_origin;


%%
alpha = 0; % AoA in degrees

Re = 2e4;
Uinf = Re .* 1.56e-5 / c;

alpha = deg2rad(alpha);
Gamma = -4 .* pi .* R .* Uinf .* sin(alpha);

%% Flow Field

n = 200
bound = 2*R;
[eta, xi] = meshgrid(linspace(-bound,bound,n), linspace(bound,-bound,n));

zeta = eta + 1i .* xi;


w = flow_field_gamma(alpha, Uinf, a, 0, zeta, Zeta_origin);

psi = imag(w);

% figure
% h = pcolor(eta,xi,psi);
% colormap hot;
% set(h, 'EdgeColor', 'none');

psi_p = sort(psi(psi>0));

nlevels = 4;

levels = imag(flow_field_gamma(alpha, Uinf, a, 0, 1i.*linspace(0,bound,nlevels), Zeta_origin));

% levels = logspace(log(psi_p(1,1)),log(psi_p(end,1)),nlevels);
levels = sort(cat(2,-levels(2:end),levels))


figure
contour(eta,xi,psi,levels)
% hold on
% plot(Zeta_circ,'b')
axis square

figure 
mesh(eta,xi,psi)



%% Plot

figure,
subplot(1,2,1), hold all
plot(complex(0),'k*') % this is (0,0)
plot(complex(Zeta_origin),'bo') % this is the center of the circle 
plot(Zeta_circ,'b')
% hold on
% plot(complex(reverse_joukowski(naca_mapping(c, tc),a)),'r')
axis equal tight, axis([-1.5*R,1.5*R,-1.5*R,1.5*R]) 
grid on, grid minor
xlabel('\eta'), ylabel('\xi')
%
subplot(1,2,2), hold all
plot(joukowski_mapping(Zeta_circ,a),'b')
% hold on
% plot(naca_mapping(c, tc),'r')
view(2), axis equal tight, axis([-0.75*c,0.75*c,-0.75*c,0.75*c]) 
grid on, grid minor
xlabel('x'), ylabel('y')

% subplot(1,3,3), hold all
% plot(naca_mapping(c, tc),'m')
% view(2), axis equal tight, axis([-0.75*c,0.75*c,-0.75*c,0.75*c]) 
% grid on, grid minor
% xlabel('x'), ylabel('y') 
clc;
clear;
close all;

set(0,'DefaultLineLineWidth',1.5);
%% Airfoil Characteristics

c = 1; % Chord length
tc = 12./100; % Thickness to Chord Ratio

%% Joukowski circle 

a = (c./4) * ( 1 + 2.*tc.*4/(3.*sqrt(3)) )  .* (1 + tc.*4/(3.*sqrt(3)) )^(-2);

b = a .* tc .* 4/(3.*sqrt(3));

R = a + b;
Zeta_origin = -b;
Zeta_circ = R*exp(1i*linspace(0,2*pi,200)) + Zeta_origin;


%% Flow Field
alpha = 0; % AoA in degrees
alpha = deg2rad(alpha);

Re = 2e4;
uinf = Re .* 1.56e-5 / c;
cl = 2*pi*sin(alpha)*R/a

Gamma = -4 .* pi .* R .* uinf .* sin(alpha); % alpha


n = 400;
bound = 3*R;
[eta, xi] = meshgrid(linspace(-bound,bound,n), linspace(bound,-bound,n));

zeta = eta + 1i .* xi;


w = flow_field_gamma(alpha, uinf, R, Gamma, zeta, Zeta_origin);

psi = imag(w);

nlevels = 6;

levels = imag(flow_field_gamma(alpha, uinf, R, Gamma, -bound + 1i.*linspace(0,bound,nlevels), Zeta_origin));

if ~ismember(0,levels)
    levels = cat(2,0,levels);
end

levels = sort(cat(2,-levels(2:end),levels));

figure
fill(real(Zeta_circ),imag(Zeta_circ),'b','LineStyle','none')
hold on
contour(eta,xi,psi,levels,'k-')
hold on
plot(complex(Zeta_origin),'kx')
axis square
grid on, grid minor
savefig('FlowContourJouk.fig')

%%
narrows = nlevels*3;

[eta2, xi2] = meshgrid(linspace(-bound,bound,narrows), linspace(bound,-bound,narrows));

zeta2 = eta2 + 1i.*xi2;

zeta2(abs(zeta2-Zeta_origin)<=R) = [];

eta2 = real(zeta2);
xi2 = imag(zeta2);

U_zeta = velocity_field_gamma(alpha, uinf, R, Gamma, zeta2, Zeta_origin);
u_zeta = real(U_zeta);
v_zeta = -imag(U_zeta);

figure
fill(real(Zeta_circ),imag(Zeta_circ),'b','LineStyle','none')
hold on
plot(complex(Zeta_origin),'kx')
hold on
quiver(eta2,xi2,u_zeta,v_zeta)
axis square
grid on, grid minor
savefig('FlowArrowsJouk.fig')

%%
[x, y] = meshgrid(linspace(-1,1,n), linspace(1,-1,n));
zed = x + 1i.*y;

[zetap, zetam] = reverse_joukowski(zed,a);

wp = flow_field_gamma(alpha, uinf, a, 0, zetap , Zeta_origin);
wm = flow_field_gamma(alpha, uinf, a, 0, zetam , Zeta_origin);

psip = imag(wp);
psim = imag(wm);

% figure 
% mesh(x,y,psip)
% hold on
% mesh(x,y,psim)
% axis square
% 
% figure 
% mesh(x,y,psip+psim)
% axis square

figure
plot(complex(zetap),'rx')
hold on
plot(complex(zetam),'bx')

% figure
% h = pcolor(x,y,t);
% axis square
% set(h, 'EdgeColor', 'none');
% 
% figure 
% mesh(x,y,t)
% axis square


levels2 = 0;
figure
plot(joukowski_mapping(Zeta_circ,a),'b')
hold on
contour(x,y,imag(joukowski_mapping(w,a)),'k-')


%% 

figure,
subplot(1,2,1), hold all
plot(complex(0),'k*') % this is (0,0)
plot(complex(Zeta_origin),'bo') % this is the center of the circle 
plot(complex(Zeta_circ),'b')
axis equal tight, axis([-1.5*R,1.5*R,-1.5*R,1.5*R]) 
grid on, grid minor
xlabel('\eta'), ylabel('\xi')


subplot(1,2,2), hold all
plot(joukowski_mapping(Zeta_circ,a),'b')
view(2), axis equal tight, axis([-0.75*c,0.75*c,-0.75*c,0.75*c]) 
grid on, grid minor
xlabel('x'), ylabel('y')

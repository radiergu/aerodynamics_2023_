clc;
clear;
close all;
%%
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0,'defaulttextinterpreter','latex');  
set(0,'defaultAxesFontSize',20)
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


figure,
subplot(1,2,1), hold all
plot(complex(Zeta_origin),'bx') % this is the center of the circle 
plot(complex(Zeta_circ),'b')
axis equal tight, axis([-1.5*R,1.5*R,-1.5*R,1.5*R]) 
grid on, grid minor
xlabel('$\eta$'), ylabel('$\xi$')

subplot(1,2,2), hold all
plot(joukowski_mapping(Zeta_circ,a),'b')
view(2), axis equal tight, axis([-0.75*c,0.75*c,-0.75*c,0.75*c]) 
grid on, grid minor
xlabel('$x$'), ylabel('$y$')



%% Flow Field

alpha = 7; % AoA in degrees
alpha = deg2rad(alpha);

backward = true;

if backward
    alpha = pi - alpha; % backward facing
end

Re = 2e4;
uinf = Re .* 1.56e-5 / c;

Gamma = -4 .* pi .* R .* uinf .* sin(alpha); % alpha

n = 500;

nacaplotflow(a,R,Zeta_origin,c,alpha,Re,n,true)
nacaplotflow(a,R,Zeta_origin,c,alpha,Re,n,false)

%%



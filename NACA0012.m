clc;
clear;
close all;
%%
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
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

Zeta0 = -b;

Zeta_circ = R*exp(1i*linspace(0,2*pi,200)) + Zeta0;


figure,
subplot(1,2,1), hold all
plot(complex(Zeta0),'bx') % this is the center of the circle 
plot(complex(Zeta_circ),'b')
axis equal tight, axis([-1.5*R,1.5*R,-1.5*R,1.5*R]) 
grid on, grid minor
xlabel('$\eta$'), ylabel('$\xi$')

subplot(1,2,2), hold all
plot(joukowski_mapping(Zeta_circ,1.*a),'b')
view(2), axis equal tight, axis([-0.75*c,0.75*c,-0.75*c,0.75*c]) 
grid on, grid minor
xlabel('$x$'), ylabel('$y$')



%% Flow Field

alpha = 7; % AoA in degrees


Re = 5*2e4;
uinf = Re .* 1.56e-5 / c;

n = 500;

nacaplotflow(a,b,c,tc,alpha,Re,n,false)
nacaplotflow(a,b,c,tc,alpha,Re,n,true)
%%

Zeta0_f = -b;
Zeta0_b = b;

Zeta_stag_f = Zeta0_f + R;
Zeta_stag_b = Zeta0_b + R;

load('naca_alpha-cl.mat');
[~, tmp1] = max(cl_nforward);
[~, tmp2] = max(cl_nreverse);
alphamax = min(alpha_nforward(tmp1), alpha_nreverse(tmp2) );
alphas = linspace(0,alphamax,50);
figure
hold on
plot(alpha_nforward, cl_nforward,'ob','DisplayName','Exp. NACA0012')
plot(alpha_nreverse, cl_nreverse,'^k','DisplayName','Exp. Reversed')
plot(alphas,lift_coeff_potential_flow(c,uinf,Circ_gamma(deg2rad(alphas),R,uinf,Zeta_stag_f,Zeta0_f)),'-r','DisplayName','potential flow theory')
hold off
grid on, grid minor
xlabel('Angle of Attack $\alpha[^\circ]$'), ylabel('Coefficient of lift $C_l$')
legend('Location','southeast')
savefig("LiftCoeff_potentialFlow.fig")
saveas(gcf,"LiftCoeff_potentialFlow",'epsc')
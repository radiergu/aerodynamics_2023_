clc;
clear;
close all;
%% Airfoil Characteristics

c = 1; % Chord length
tc = 12./100; % Thickness to Chord Ratio

%% Joukowski circle 

a = (c./4) * ( 1 + 2.*tc.*4/(3.*sqrt(3)) )  .* (1 + tc.*4/(3.*sqrt(3)) )^(-2)

b = a .* tc .* 4/(3.*sqrt(3))

R = a + b;
Zeta_origin = b;
Zeta_circ = R*exp(1i*linspace(0,2*pi,100)) + Zeta_origin;

%% Plot

figure, 
subplot(1,2,1), hold all
plot(complex(0),'k*') % this is (0,0)
plot(complex(Zeta_origin),'bo') % this is the center of the circle 
plot(Zeta_circ,'b')
hold on
plot(complex(reverse_joukowski(naca_mapping(c, tc),a)),'r')
axis equal tight, axis([-1.5*R,1.5*R,-1.5*R,1.5*R]) 
grid on, grid minor
xlabel('\eta'), ylabel('\xi')
%
subplot(1,2,2), hold all
plot(joukowski_mapping(Zeta_circ,a),'b')
hold on
plot(naca_mapping(c, tc),'r')
view(2), axis equal tight, axis([-0.75*c,0.75*c,-0.75*c,0.75*c]) 
grid on, grid minor
xlabel('x'), ylabel('y')

% subplot(1,3,3), hold all
% plot(naca_mapping(c, tc),'m')
% view(2), axis equal tight, axis([-0.75*c,0.75*c,-0.75*c,0.75*c]) 
% grid on, grid minor
% xlabel('x'), ylabel('y') 
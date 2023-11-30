clc;
clear;
close all;

%% generate airfoil

c = 1; % Chord length
tc = 9./100; % Thickness to Chord Rati

z = camber_plate(c, tc, 'yes');

%% calculation of the lift coefficient

alpha = linspace(-3, 15, 20);

% for the camberline we considered we have 

A0 = deg2rad(alpha);
A1(1:length(A0)) = 4*tc;

cl = pi.*(2 * A0 + A1); 

figure
plot(alpha, cl, 'b')

%% checking severall camb


tc_ = floor(linspace(0,9,3))./100;

figure; hold on; grid on;
legend = [];

A0 = deg2rad(alpha);

for k = 1:3
    
    A1(1:length(A0)) = 4*tc_(k);

    cl = pi.*(2 * A0 + A1);

    plot(alpha, cl)
    legend(k) = "tc = " + string(tc);
    
end

legend(legend);



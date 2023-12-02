clc;
clear;
close all;

%% Open some datas

[numbers, ~] = xlsread("Extracted data.xlsx", 'Feuille 1'); %#ok<XLSRD>
if ~isempty(numbers)
    data =  numbers;
end

x_3per = rmmissing(data(:,6));
y_3per = rmmissing(data(:,7));

x_9per = rmmissing(data(:,8));
y_9per = rmmissing(data(:,9));

% list 
alpha_paper = {x_3per, x_9per};
cl_paper    = {y_3per, y_9per};

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

tc_ = [3, 9]./100;

figure(2); hold on; grid on;
legend_ = [""];

A0 = deg2rad(alpha);

l = 1;

for k = 1:2
    
    A1(1:length(A0)) = 4*tc_(k);

    cl = pi.*(2 * A0 + A1);

    plot(alpha, cl)
    plot(alpha_paper{k}, cl_paper{k})

    legend_(l)     = "tc = " + string(tc);
    legend_(l + 1) = "Paper, tc = " + string(tc);
    l = l + 2;
end
    
    figure(2);
    legend(legend_);



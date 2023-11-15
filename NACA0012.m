%% Airfoil Characteristics

c = 1; % Chord length
tc = 12./100; % Thickness to Chord Ratio

a = (c./4) * ( 1 + 2.*tc.*4/(3.*sqrt(3)) )  .* (1 + tc.*4/(3.*sqrt(3)) )^(-2)

b = a .* tc .* 4/(3.*sqrt(3))


%% Joukowski circle 

R = a + b


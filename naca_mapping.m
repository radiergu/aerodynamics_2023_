function z = naca_mapping(c, tc)
x = linspace(0,c,200);

yplus = 5.*tc.*(0.2969.*sqrt(x) - 0.1260.*x - 0.3516.*x.^2 + 0.2843.*x.^3 - 0.1036.*x.^4 );
yminus = - yplus;

z = cat(2, -(x-c./2) + 1i.*yplus, flip(-(x-c./2) + 1i.*yminus) );
end
function nacaplotflow(a,R,Zeta0,c,alpha,Re,n,backward)

Zeta_circ = R*exp(1i*linspace(0,2*pi,200)) + Zeta0;

if backward
    alpha = pi - alpha; % backward facing
end

uinf = Re .* 1.56e-5 / c;

Gamma = -4 .* pi .* R .* uinf .* sin(alpha); % alpha

bound = 3*R;


[eta, xi] = meshgrid(linspace(-bound,bound,n), linspace(bound,-bound,n));
zeta = eta + 1i .* xi;

cond = abs(zeta-Zeta0)<=R;
zeta(cond) = (R).*exp(-1i.* linspace(0,2*pi,length(zeta(cond))) ) + Zeta0;

eta = real(zeta);
xi = imag(zeta);

x = real(joukowski_mapping(zeta,a));
y = imag(joukowski_mapping(zeta,a));

w = flow_field_gamma(alpha, uinf, R, Gamma, zeta, Zeta0);

psi = imag(w);

U_zeta = velocity_field_gamma(alpha, uinf, R, Gamma, zeta, Zeta0);

% cond2 = abs((x+1i.*y) )<=1/(1.5*100);
% U_zeta(cond) = 0;
U_real = U_zeta ./(1 - a.^2./(zeta.^2) );
Umod = abs(U_real);
% U_mod(cond2) = 0;

% cond = abs(Umod)> mean(Umod,"all") + 2*std(Umod(:));
% Umod(cond) = 0; %mean(Umod,"all") + 2*std(Umod(:));



nlevels = 6;

levels = imag(flow_field_gamma(alpha, uinf, R, Gamma, (-10*R + 1i.*linspace(0,1.5*R,nlevels)).*exp(1i.*alpha) , Zeta0));

if ~ismember(0,levels)
    levels = cat(2,0,levels);
end

levels = sort(cat(2,-levels(2:end),levels));

figure
hold on
colormap hot;
h1 = pcolor(eta,xi,abs(U_zeta));
contour(eta,xi,psi,levels,'k-')
fill(real(Zeta_circ),imag(Zeta_circ),'b','LineStyle','-')
hold off
axis square
grid on, grid minor
set(h1, 'EdgeColor', 'none');
colorbar
name = "FlowContourJouk";
if backward
    name = name + "_backward";
else
    name = name + "_forward";
end
savefig(name + '.fig')



figure
hold on
colormap hot;
%h2 = pcolor(x,y,abs(U_real));
h2 = pcolor(x,y,Umod);
contour(x,y,psi,levels,'k-')
fill(real(joukowski_mapping(Zeta_circ,a)),imag(joukowski_mapping(Zeta_circ,a)),'b','LineStyle','-')
hold off
axis square
grid on, grid minor
set(h2, 'EdgeColor', 'none');
colorbar
name2 = "FlowContourFoil";
if backward
    name2 = name2 + "_backward";
else
    name2 = name2 + "_forward";
end
savefig(name2 + '.fig')
end
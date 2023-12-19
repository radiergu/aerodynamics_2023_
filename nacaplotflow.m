function nacaplotflow(a,R,Zeta0,c,tc,alpha,Re,n,backward)

Zeta_circ = R*exp(1i*linspace(0,2*pi,200)) + Zeta0;

if backward
    % alpha = pi - alpha; % backward facing
end

uinf = Re .* 1.56e-5 / c;

Gamma = -4 .* pi .* R .* uinf .* sin(alpha);

if backward
    Gamma = 1i.*2.*pi.*uinf .* (a.*exp(-1i.*alpha) ...
        - R.^2.*exp(1i.*alpha)./a );
    Gamma = -4 .* pi .* R .* uinf .* sin(alpha);
end



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
U_real = U_zeta ./(1 - a.^2./(zeta.^2) );
Umod = abs(U_real);



nlevels = 6;

levels = imag(flow_field_gamma(alpha, uinf, R, Gamma, (-10*R + 1i.*linspace(0,1.5*R,nlevels)).*exp(1i.*alpha) , Zeta0));

if ~ismember(0,levels)
    levels = cat(2,0,levels);
end

levels = sort(cat(2,-levels(2:end),levels));

figure
hold on
colormap turbo;
h1 = pcolor(eta,xi,abs(U_zeta));
contour(eta,xi,psi,levels,'k-')
fill(real(Zeta_circ),imag(Zeta_circ),'w','LineStyle','-')
hold off
lims = [-2 2].*R+Zeta0;
xlim(lims)
ylim(lims)
axis square
grid on, grid minor
set(h1, 'EdgeColor', 'none');
colorbar
clim([min(abs(U_zeta),[],"all"), max(abs(U_zeta),[],"all")]);
name = "FlowContourJouk";
if backward
    name = name + "_backward";
else
    name = name + "_forward";
end
savefig(name + '.fig')



figure
hold on
colormap turbo;
%h2 = pcolor(x,y,abs(U_real));
h2 = pcolor(x,y,Umod);
contour(x,y,psi,levels,'k-')
fill(real(joukowski_mapping(Zeta_circ,a)),imag(joukowski_mapping(Zeta_circ,a)),'w','LineStyle','-')
hold off
xlim([-0.75 0.75].*c)
ylim([-4 4].*tc*c)

grid on, grid minor
set(h2, 'EdgeColor', 'none');
colorbar
clim([min(Umod,[],"all"), max(Umod,[],"all")]);
name2 = "FlowContourFoil";
if backward
    % set(gca,'xdir','reverse')
    name2 = name2 + "_backward";
else
    name2 = name2 + "_forward";
end
savefig(name2 + ".fig")
end
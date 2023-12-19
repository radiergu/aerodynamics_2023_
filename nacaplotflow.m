function nacaplotflow(a,b,c,tc,alpha,Re,n,backward)
alpha2 = deg2rad(alpha);
R = a + b;

Zeta0 = -b;
if backward
    Zeta0 = b;
end

Zeta_circ = R*exp(1i*linspace(0,2*pi,200)) + Zeta0;

Zeta_stag = Zeta0 + R;

if backward
    Zeta_stag = Zeta0 - R;
end



uinf = Re .* 1.56e-5 / c;

Gamma = Circ_gamma(alpha2, R, uinf, Zeta_stag, Zeta0);




bound = R*3;


[eta, xi] = meshgrid(linspace(-bound,bound,n), linspace(bound,-bound,n));
zeta = eta + 1i .* xi;

cond = abs(zeta-Zeta0)<=R;
zeta(cond) = (R).*exp(-1i.* linspace(0,2*pi,length(zeta(cond))) ) + Zeta0;

eta = real(zeta);
xi = imag(zeta);

x = real(joukowski_mapping(zeta,a));
y = imag(joukowski_mapping(zeta,a));

w = flow_field_gamma(alpha2, uinf, R, Gamma, zeta, Zeta0);

psi = imag(w);

U_zeta = velocity_field_gamma(alpha2, uinf, R, Gamma, zeta, Zeta0);
U_real = U_zeta ./(1 - a.^2./(zeta.^2) );
Umod = abs(U_real);



nlevels = 4;

levels = imag(flow_field_gamma(alpha2, uinf, R, Gamma, (-10*R + 1i.*linspace(0,1.5*R,nlevels)).*exp(1i.*alpha2) , Zeta0));

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
name = "FlowContourJouk_alpha="+string(alpha);
if backward
    name = name + "_backward";
else
    name = name + "_forward";
end
savefig(name + '.fig')
saveas(gcf,name,'epsc')


figure
hold on
colormap turbo;
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
name2 = "FlowContourFoil_alpha="+string(alpha);
if backward
    name2 = name2 + "_backward";
else
    name2 = name2 + "_forward";
end
savefig(name2 + ".fig")
saveas(gcf,name2,'epsc')



end
function w = flow_field_gamma(alpha, Uinf, a, Gamma, zeta)

w = Uinf .*( zeta .* exp(-1i.*alpha) + a.^2./zeta .* exp(1i.*alpha) ) - 1i .* Gamma ./(2.*pi) .* log(zeta);

end
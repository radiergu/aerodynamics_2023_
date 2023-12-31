function U = velocity_field_gamma(alpha, Uinf, R, Gamma, zeta, zeta0)

U = Uinf .* exp(-1i.*alpha) ...
    - Uinf .* R.^2 .* exp(1i.*alpha) .* (zeta - zeta0).^(-2) ...
    - 1i.*Gamma./(2.*pi .* (zeta-zeta0));
end
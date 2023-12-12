function w = flow_field_gamma(alpha, Uinf, R, Gamma, zeta, zeta0)

w = Uinf .*( (zeta - zeta0) .* exp(-1i.*alpha) ...
    + R.^2./(zeta - zeta0) .* exp(1i.*alpha) ) ...
    - 1i .* Gamma ./(2.*pi) .* log(zeta - zeta0);

end
function cl = lift_coeff_potential_flow(c, uinf, Gamma)
    cl = -2.*Gamma ./ (uinf .* c);
end
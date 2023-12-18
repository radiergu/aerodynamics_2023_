function cl = lift_coeff_potential_flow(alpha, tc)
    cl = 2*pi*sin(alpha).*(1 + tc.*4/(3.*sqrt(3)) );
end
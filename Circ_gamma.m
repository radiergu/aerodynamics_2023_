function Gamma = Circ_gamma(alpha, R, uinf, Zeta_stag, Zeta0)
    Gamma = -1i.*2.*pi.*uinf.*(exp(-1i.*alpha).*(Zeta_stag - Zeta0) - R^2.*exp(1i.*alpha)./(Zeta_stag - Zeta0)  );
end
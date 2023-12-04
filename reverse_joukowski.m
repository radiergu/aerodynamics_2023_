function zeta = reverse_joukowski(z, a)
zeta = z./2 + ( (z./2).^2 - a.^2 ).^(1/2);
%zeta = cat(2, zeta, flip(-zeta) );
end
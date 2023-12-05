function [zetap, zetam] = reverse_joukowski(z, a)
zetap = z./2 + ( (z.^2)./4 - a.^2 ).^(1/2);
zetam = z./2 - ( (z.^2)./4 - a.^2 ).^(1/2);
end
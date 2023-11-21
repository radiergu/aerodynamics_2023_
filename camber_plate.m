function [z] = camber_plate(c, tc, do_plot)
% Generate a camber line with cord length c and t the percentage of c 
% expected for the height at c/2.
%
% here xc = x/c
%
% dz is the vector of complex number  of the form x + i*dyc

xc = linspace(0,1,200); 

yc  = 4 * tc * xc .* (1 - xc);
dycdx = 4 * tc .* (1 - 2 * xc);     

z  = complex(xc * c + 1i * yc);

if do_plot == "yes"
    
    figure; hold on; grid on; xlim([0 c]); ylim([-c c]); 
    plot(z, 'b', 'LineWidth', 4)
    yline(0,'k', 'LineStyle','--');
  
end 
end





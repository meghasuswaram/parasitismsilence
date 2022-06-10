function [ Rpar ] = repropar( x, zbar, rmin, rmax, beta )
diff_in_trait = (x - zbar);
denom = (1 + exp(-(beta.*diff_in_trait)));
second_term = (rmax - rmin)./denom;
Rpar = rmin + second_term;


end




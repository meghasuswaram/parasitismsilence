function [ R ] = repro( x, zbar, rmin, rmax, alphap , nz)
alpha = (alphap/sum(nz));
diff_in_trait = (x - zbar);
denom = (1 + exp(-(alpha.*diff_in_trait)));
second_term = (rmax - rmin)./denom;
R = rmin + second_term;

 % R = rmin + ((rmax - rmin)./(1 + exp(-((alphap./nz).*(x-zbar)))));
% zbar
% alphap
% sum(nz)
% rmax
% rmin

end




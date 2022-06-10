function [ Probescpar1 ] = probescpar( x, nz, pt, fmax)
denomparesc = fmax + (sum(nz) .* exp(-(x/sum(nz))));
numerpar = exp(-(x/sum(nz))).*sum(pt)*fmax;

Probescpar1 = exp(-(denomparesc./numerpar));

end




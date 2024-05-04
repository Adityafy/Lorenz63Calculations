function mat = lor63jacobian(p,Xt)
    sigma = p.sigma;
    r = p.r;
    b = p.b;
    mat = [-sigma sigma 0; ...
        (r-Xt(3)) -1 -Xt(1); ...
        Xt(2) Xt(1) -b];
end

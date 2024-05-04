function L = linLorenz63(p)
sig = p.sigma;
r = p.r;
b = p.b;

L = [-sig, sig, 0;
        r, -1, 0;
        0, 0, -b];
end
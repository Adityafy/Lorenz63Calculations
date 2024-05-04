function mat = rosslerjacobian(p,Xt)
    a = p.a;
    b = p.b;
    c = p.c;
    mat = [0 -1 -1; ...
        1 a 0; ...
        Xt(3) 0 Xt(1)-c];
end

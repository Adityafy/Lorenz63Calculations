function mat = lor63jacobian(para,Xt,vt)
    sig = para(1);
    r = para(2);
    b = para(3);
    mat = [-sig sig 0; ...
        (r-Xt(3)) -1 -Xt(1); ...
        Xt(2) Xt(1) -b] * vt;
end

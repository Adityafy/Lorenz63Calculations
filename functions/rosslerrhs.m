function vec = rosslerrhs(p,X)
    a = p.a;
    b = p.b;
    c = p.c;
    x = -X(2,:)-X(3,:);
    y = X(1,:)+a*X(2,:);
    z = b+X(3,:).*(X(1,:)-c);
    vec = [x; y; z];
end
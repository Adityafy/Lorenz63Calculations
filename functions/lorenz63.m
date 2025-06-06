function vec = lorenz63(p,X)
    sig = p.sig;
    r = p.r;
    b = p.b;
    x = sig*(X(2,:)-X(1,:));
    y = r*X(1,:)-X(2,:)-X(1,:)*X(3,:);
    z = X(1,:)*X(2,:)-b*X(3,:);
    vec = [x; y; z];
end
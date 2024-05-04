function vec = lorenz63rhs(p,X)
    sigma = p.sigma;
    r = p.r;
    b = p.b;
    x = sigma*(X(2,:)-X(1,:));
    y = r*X(1,:)-X(2,:)-X(1,:)*X(3,:);
    z = X(1,:)*X(2,:)-b*X(3,:);
    vec = [x; y; z];
end
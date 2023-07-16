function vec = lorenz63(para,X)
    sig = para(1);
    r = para(2);
    b = para(3);
    x = sig*(X(2,:)-X(1,:));
    y = r*X(1,:)-X(2,:)-X(1,:)*X(3,:);
    z = X(1,:)*X(2,:)-b*X(3,:);
    vec = [x; y; z];
end
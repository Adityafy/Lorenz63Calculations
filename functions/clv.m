function [wf,cle] = clv(p,dynamicsMat,gsVec,cmat,evolFunc,jacobianFunc)
% [wf,cle] = clv(p,dynamicsMat,gsVec,cmat,evolFunc,jacobianFunc)

[M,total_time] = size(dynamicsMat);
wi = zeros(M,M, total_time); % w (clv) initial
% wmag(:,1) = zeros(Nx,1); % w magnitude
cle(:,1) = zeros(M,1);

fprintf("\nCalculating CLVs...\n")
for t = 1:total_time
    for j = 1:M
        for k = 1:j
            wi(:,j,t) = wi(:,j,t) + cmat(k,j,t).*gsVec(:,k,t);
        end
    end
    if rem(t,p.wnorm) == 0
        for j = 1:M
            wi(:,j,t) = wi(:,j,t)./norm(wi(:,j,t));
        end
    end
    if rem(t, total_time/10) == 0
        fprintf("This is time step %g and ", t);
        toc
    end
end

wf = zeros(M, M, total_time); % w final
for t = 1:total_time
    wf(:,:,t) = evolFunc(p,dynamicsMat(:,t),wi(:,:,t),jacobianFunc);
    for j = 1:M
        %wmag(j,t) = norm(wf(:,j,t)); % being normalized at every time step
        %(not explicitely normalized)
        wmag = norm(wf(:,j,t));
        %cle(j,t) = log(abs(wmag(j,t)));
        cle(j,t) = (1/p.twnorm)*log(abs(wmag));
        wf(:,j,t) = wf(:,j,t)./wmag;
    end
end
toc
end
function u = stdynamics(p,u,solverfunc,dynrhsfunc,linpartfunc,nonlinpartfunc)
% 
% calculates the space time dynamics

u = solverfunc(p,u,dynrhsfunc,linpartfunc,nonlinpartfunc);
end
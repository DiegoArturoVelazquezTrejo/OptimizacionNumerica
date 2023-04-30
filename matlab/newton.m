function [pk] = newton(xk, gradkxk, hessiana)
% Resolución del problema a través del método de newton. 
hk = hessiana(xk); 
pk = linsolve(hk, -gradkxk); 
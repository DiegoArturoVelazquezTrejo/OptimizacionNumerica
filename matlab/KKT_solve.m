function [x, lambda] = KKT_solve(G, c, A, b, x0, metodo)
% Método que resuelve un problema cuadrático utilizando KKT 
n = size(A, 2); 
p = size(A, 1); 
% Construímos la matriz de KKT 
K = [G, A.'; A, zeros(p)]; 
h = A*x0 - b; 
g = c + G * x; 
% Construímos el vector con el que resolveremos el sistema de ecuaciones 
m = [g; h]; 

if(metodo == "linsolve")
    linsolve(K, m); 
end 
if(metodo == "cholesky")
    
end 
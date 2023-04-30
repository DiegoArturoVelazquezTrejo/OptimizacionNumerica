function [NORMA, iters] = GradienteConjugadoVarianteB(A, b, x0, EPS)

% Inicialización de los parámetros
rk = A*x0 -b; 
pk = -rk; 
xk = x0; 
% Vector con todas las x_stars 
%XKS = []; 
% Iteraciones 
i = 1; 
%XKS(:,1) = xk;

% Arreglo para guardar la norma de xk 
NORMA = []; 

while(norm(rk) > EPS)

    % Obtenemos el tamaño de paso 
    alphak = (rk'*rk)/(pk'*A*pk);  
    xk = xk + alphak * pk; 
    rk_new = rk + alphak * A * pk; 
    % Cálculo de la Bk 
    Bk = (rk_new'*rk_new)/(rk'*rk);  
    rk = rk_new; 
    % Cálculo del paso 
    pk = -rk + Bk*pk;
    % Vamos guardando el vector 
    i = i + 1;
    %XKS(:,i) = xk; 
    NORMA(i) = norm(rk); 
end 
iters = i; 
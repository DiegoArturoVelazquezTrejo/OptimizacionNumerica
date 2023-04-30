function [NORMA, iters] = GradienteConjugadoPrecondicionado(A, b, x0, M, EPS)

% Inicialización de los parámetros
rk = A*x0 -b; 

% Vamos a resolver el siguiente sistema de ecuaciones: M y0 = r0 
yk = solve(M, rk);

pk = -yk;  
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
    alphak = (rk'*yk)/(pk'*A*pk); 
    xk = xk + alphak * pk; 
    rk_new = rk + alphak * A * pk;
    % Resolvemos M yk+1 = rk+1
    yk_new = solve(M, rk_new); 
    % Cálculo de la Bk 
    Bk = (rk_new'*yk_new)/(rk'*yk); 
    % Cálculo del paso 
    pk = -yk_new + Bk*pk;

    % Vamos guardando el vector 
    i = i + 1;
    %XKS(:,i) = xk; 
    NORMA(i) = norm(rk); 
    rk = rk_new;
    yk = yk_new;
end 
iters = i; 
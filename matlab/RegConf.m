function [xk, k, XSKS, FXKS] = RegConf(f, df, x0, metodo, metodoBK, delta_hat, etha, EPS, MAXITER, B)
% Método de regiones por confianza para la optimización de una función   
xk = x0; 
Bk = eye(length(xk)); 
delta_k = 0.1; 
gk = df(x0); 
norma = norm(gk);
% Listas con los estados 
XSKS = [];
FXKS = []; 
k = 1; 
while(norma > EPS && k < MAXITER && delta_k > EPS)
    XSKS(:,k) = xk; 
    FXKS(k) = f(xk);
    % Defninimos mk 
    mk = @(p) f(xk) + gk'*p + 0.5*p'*Bk*p; 
    pk = calcula_pk(delta_k, gk, Bk, metodo); 
    rk = (f(xk) - f(xk + pk))/(mk(zeros(size(xk))) - mk(pk));

    if(rk < 1/4)
        delta_k = (1/4)* delta_k; 
    else
        if(rk > 3/4 && abs(norm(pk) - delta_k) < EPS)
            delta_k = min(2*delta_k, delta_hat); 
        else 
            delta_k = delta_k + 0; % Esto se puede omitir 
        end
    end 
    if(rk > etha)
        xk_new = xk + pk; 
        % Recalculamos Bk 
        gradxk1 = df(xk_new); 
        gradxk = df(xk);  
        if(metodo == "hessiana")
            Bk = B(xk); 
        end
        %Bk = actualizaBK(metodoBK, Bk, xk_new, xk, gradxk1, gradxk); 
    else
        xk_new = xk + 0;  % Esto también se puede omitir 
    end 
    k = k + 1; 
    xk = xk_new; 
    gk = df(xk_new); 
    norma = norm(gk); 
    circplot(xk, delta_k, "r--");
    hold on; 

end 
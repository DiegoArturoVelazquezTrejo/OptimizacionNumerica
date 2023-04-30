function [xk, k, XKS, FXKS] = BusquedaLinea(x0, f, grad, hessiana, MAXITER, EPS, metodo)
% Método de búsqueda en línea con la optimización que se puede escoger el
% metodo de optimización para hallar la pk. 
k = 0; 
xk = x0; 
BK = eye(length(xk)); 
norma = norm(grad(xk));
XKS = [];
FXKS = []; 
while(norma > EPS && k < MAXITER) 
        gradxk = grad(xk); 
        if(metodo == "GD")
            pk = GD(gradxk); 
        end
        if(metodo == "newton")
            pk = newton(xk, gradxk, hessiana); 
        end
        if(metodo == "SR1" || metodo == "BFGS")
            pk = SR1_BFGS(gradxk, BK); 
        end  
        ak = bisect(f, xk, pk);
        xk_new = xk + ak * pk;
        k = k + 1; 
        gradk1 = grad(xk_new); 
        BK = actualizaBK(metodo, BK, xk_new, xk, gradk1, gradxk); 
        norma = norm(gradk1); 
        XKS(:,k) = xk; 
        FXKS(k) = f(xk);
        xk = xk_new; 
end
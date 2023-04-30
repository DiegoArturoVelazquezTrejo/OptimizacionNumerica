function [pk] = calcula_pk(delta_k, gk, Bk, metodo)

gkBkgk = gk'*Bk*gk; 
% Método que actualiza la pk dependiendo el método 
if(metodo == "cauchy")
    norma_gk = norm(gk); 
    if(gkBkgk < 0)
        tauk = 1; 
    else 
        tauk = min(norma_gk^3/(delta_k * gkBkgk), 1); 
    end 
    ps = (-delta_k / norma_gk) * gk; 
    pk = tauk*ps; 
end 
if(metodo == "dogleg") 
    pb = linsolve(Bk, -gk);
    pu = -(gk'*gk/gkBkgk) * gk; 
    a = norm(pb - pu)^2;
    b = 2 * pu'*(pb - pu);
    c = norm(pu)^2 - delta_k^2;
    tauk = -1;
    det = b^2 - 4 * a * c;
    if (a > 1e-4 && det > 0)
        tauk = max((-b - sqrt(det))/(2*a), (-b + sqrt(det))/(2*a));
        tauk = tauk + 1;
    end
    if(tauk >= 0 && tauk <= 1)
        pk = tauk * pu; 
    elseif(tauk > 1 && tauk <= 2)
        pk = pu + (tauk - 1)* (pb - pu); 
    else
        pk = -(delta_k / norm(gk)) * gk;
    end 
end
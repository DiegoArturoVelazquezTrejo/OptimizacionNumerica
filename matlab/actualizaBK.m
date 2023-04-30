function [BK] = actualizaBK(metodo, BK, xk1, xk, gradxk1, gradxk)
% Función que modifica la matriz BK de acuerdo al método (SR1 o BFGS)
sk = xk1 - xk; 
yk = gradxk1 - gradxk; 

if(metodo == "SR1")
    AK = BK * sk; 
    YK = yk - AK; 
    num = YK * YK.'; 
    den = YK.' * sk; 
    BK = BK + num/den; 
end 

if(metodo == "BFGS")
    DK =  sk.' * BK; 
    AK = BK * sk; 
    num1 = AK * DK; 
    den1 = DK * sk; 
    num2 = yk * yk.'; 
    den2 = yk.' * sk; 
    BK =  BK - num1/den1 + num2/den2; 
end 
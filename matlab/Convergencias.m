
metodos = {"newton", "SR1", "BFGS", "GD"}; 
colores = {'k--', 'r', 'g*:', 'b--'}; 
% Parámetros globales 
MAXITER = 1000; 
EPS = 1e-20; 
% Vector inicial 
x0 = [1;1;1;1;1;1;1;1;1;1]; 
 % genera matriz aleatoria de 10x1
b = randn(10, 1);

subplot(3, 2, 1); 

for k = 1:6
    % Dimensión de la matriz 
    n = 10; 
    % Generamos una matriz positiva definida de tamaño nxn 
    U = gallery('orthog', n); 
    % Construímos la matriz V 
    L = diag((1:n).^k); 
    Q = U * L * U.';
    f = @(x) x.' * Q * x; %+ b.' * x; 
    df = @(x) 2 * Q * x;  %+ b; 
    hf = @(x) 2 * Q; 
    subplot(3, 2, k); 
    for j = 1:4
        [xk, iteraciones, XSKS, FXKS] = BusquedaLinea(x0, f, df, hf, MAXITER, EPS, metodos{j}); 
        %disp(metodos{j}, iteraciones);  
        plot(1:iteraciones,log10(FXKS),colores{j}); 
        hold on 
    end 
    legend(metodos); 
    title("Exponente k= %d", k); 
end 
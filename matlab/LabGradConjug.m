
% Inicialización de los parámetros 
n = 1000
A = gallery('gcdmat', n)
b = ones(n, 1)

x0 = rand(n, 1)
% Solución por gradiente conjugado 
[xks, NORMA, iters] = GradienteConjugado(A, b, x0, 1e-6)


figure;
subplot(2,1,1); % subtrama 1
title('Norma del residuo para n=10');

% Graficamos los residuos 
plot(log10(NORMA))

% Ahora vamos a hacerlo para distintas n 
NITER = 100; 
ITERACIONES = []
for i = 1:NITER
    % Inicialización de los parámetros 
    A = gallery('lehmer', i); 
    b = ones(i, 1); 
    
    x0 = rand(i, 1); 
    % Solución por gradiente conjugado 
    [xks, NORMA, iters] = GradienteConjugado(A, b, x0, 1e-6); 
    ITERACIONES(i) = iters; 
end 

subplot(2,1,2); % subtrama 2
title('Gráfica de iteraciones');
plot(ITERACIONES)


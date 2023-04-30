% Inicialización de los parámetros 
n = 100; 

Q = gallery('orthog', n); 
L1 = diag(1:n); 
A1 = Q'*L1*Q; 

k = round(n/2); 
L2 = diag( [1:k n*ones(1, n-k)]); 
A2 = Q'*L2*Q;

L3 = diag([1:k n*ones(1, n-k) + rand(1, n-k)]);
A3 = Q'*L3*Q; 

% Condiciones iniciales
b = ones(n, 1); 
x0 = rand(n, 1); 

% Gráfica para la norma 
[NORMA1, iters1]    = GradienteConjugadoVarianteB(A1, b, x0, 1e-6); 
[NORMA2, iters2]    = GradienteConjugadoVarianteB(A2, b, x0, 1e-6); 
[NORMA3, iters3]    = GradienteConjugadoVarianteB(A3, b, x0, 1e-6); 

figure;
subplot(2,1,1); % subtrama 1

plot(log10(NORMA1)); 
hold on; 
plot(log10(NORMA2));
plot(log10(NORMA3));
legend('Uniform', 'Clustered1', 'Clustered2'); 
title('Norma del residuo para variante GC');

hold on; 

% Lo mismo pero para GC 
% Gráfica para la norma 
[NORMA1, iters1]    = GradienteConjugado(A1, b, x0, 1e-6); 
[NORMA2, iters2]    = GradienteConjugado(A2, b, x0, 1e-6); 
[NORMA3, iters3]    = GradienteConjugado(A3, b, x0, 1e-6); 

subplot(2,1,2); % subtrama 1

plot(log10(NORMA1)); 
hold on; 
plot(log10(NORMA2)); 
plot(log10(NORMA3));
legend('Uniform', 'Clustered1', 'Clustered2'); 
title('Norma del residuo para GC');

% Parámetro global 
NITER = 500; 


% ------------------------------ GRADIENTE CONJUGADO
% Iteraciones 
% Ahora vamos a hacerlo para distintas n 
ITERACIONES1 = []; 
ITERACIONES2 = []; 
ITERACIONES3 = []; 

for i = 1:NITER
    
    % Construimos las variables 
    Q = gallery('orthog', i); 
    L1 = diag(1:i); 
    A1 = Q'*L1*Q; 
    
    k = round(i/2); 
    L2 = diag( [1:k i*ones(1, i-k)]); 
    A2 = Q'*L2*Q;
    
    L3 = diag([1:k i*ones(1, i-k) + rand(1, i-k)]);
    A3 = Q'*L3*Q; 
    
    % Condiciones iniciales
    b = ones(i, 1); 
    x0 = rand(i, 1); 

    % Solución por gradiente conjugado 
    [NORMA1, iters1] = GradienteConjugado(A1, b, x0, 1e-6); 
    ITERACIONES1(i) = iters1;

    % Solución por gradiente conjugado 
    [NORMA2, iters2] = GradienteConjugado(A2, b, x0, 1e-6); 
    ITERACIONES2(i) = iters2;

    % Solución por gradiente conjugado 
    [NORMA3, iters3] = GradienteConjugado(A3, b, x0, 1e-6); 
    ITERACIONES3(i) = iters3;

end 

figure; 
plot(ITERACIONES1); 
hold on; 
plot(ITERACIONES2); 
plot(ITERACIONES3); 
legend('Uniform', 'Clustered1', 'Clustered2');
title('Iteraciones vs n : Gradiente Conjugado'); 

% Tiempo  
% Ahora vamos a hacerlo para distintas n 
TIEMPO1 = []; 
TIEMPO2 = []; 
TIEMPO3 = []; 

for i = 1:NITER
    
    % Construimos las variables 
    Q = gallery('orthog', i); 
    L1 = diag(1:i); 
    A1 = Q'*L1*Q; 
    
    k = round(i/2); 
    L2 = diag( [1:k i*ones(1, i-k)]); 
    A2 = Q'*L2*Q;
    
    L3 = diag([1:k i*ones(1, i-k) + rand(1, i-k)]);
    A3 = Q'*L3*Q; 
    
    % Condiciones iniciales
    b = ones(i, 1); 
    x0 = rand(i, 1); 
    
    % Solución por gradiente conjugado 
    tic; 
    [NORMA1, iters1] = GradienteConjugado(A1, b, x0, 1e-6); 
    TIEMPO1(i) = toc;

    % Solución por gradiente conjugado 
    tic 
    [NORMA2, iters2] = GradienteConjugado(A2, b, x0, 1e-6); 
    TIEMPO2(i) = toc;

    % Solución por gradiente conjugado 
    tic; 
    [NORMA3, iters3] = GradienteConjugado(A3, b, x0, 1e-6); 
    TIEMPO3(i) = toc;

end 

figure; 
plot(TIEMPO1); 
hold on; 
plot(TIEMPO2); 
plot(TIEMPO3); 
legend('Uniform', 'Clustered1', 'Clustered2'); 
title('TIEMPO vs n : Gradiente Conjugado'); 


% ------------------------------ GRADIENTE CONJUGADO VARIANTE
% Iteraciones 
% Ahora vamos a hacerlo para distintas n 
ITERACIONES1 = []; 
ITERACIONES2 = []; 
ITERACIONES3 = []; 

for i = 1:NITER
    
    % Construimos las variables 
    Q = gallery('orthog', i); 
    L1 = diag(1:i); 
    A1 = Q'*L1*Q; 
    
    k = round(i/2); 
    L2 = diag( [1:k i*ones(1, i-k)]); 
    A2 = Q'*L2*Q;
    
    L3 = diag([1:k i*ones(1, i-k) + rand(1, i-k)]);
    A3 = Q'*L3*Q; 
    
    % Condiciones iniciales
    b = ones(i, 1); 
    x0 = rand(i, 1); 

    % Solución por gradiente conjugado 
    [NORMA1, iters1] = GradienteConjugadoVarianteB(A1, b, x0, 1e-6); 
    ITERACIONES1(i) = iters1;

    % Solución por gradiente conjugado 
    [NORMA2, iters2] = GradienteConjugadoVarianteB(A2, b, x0, 1e-6); 
    ITERACIONES2(i) = iters2;

    % Solución por gradiente conjugado 
    [NORMA3, iters3] = GradienteConjugadoVarianteB(A3, b, x0, 1e-6); 
    ITERACIONES3(i) = iters3;

end 

figure; 
plot(ITERACIONES1); 
hold on; 
plot(ITERACIONES2); 
plot(ITERACIONES3); 
legend('Uniform', 'Clustered1', 'Clustered2');
title('Iteraciones vs n : Gradiente Conjugado VARIANTE'); 


% Tiempo  
% Ahora vamos a hacerlo para distintas n 
TIEMPO1var = []; 
TIEMPO2var = []; 
TIEMPO3var = []; 

for i = 1:NITER
    
    % Construimos las variables 
    Q = gallery('orthog', i); 
    L1 = diag(1:i); 
    A1 = Q'*L1*Q; 
    
    k = round(i/2); 
    L2 = diag( [1:k i*ones(1, i-k)]); 
    A2 = Q'*L2*Q;
    
    L3 = diag([1:k i*ones(1, i-k) + rand(1, i-k)]);
    A3 = Q'*L3*Q; 
    
    % Condiciones iniciales
    b = ones(i, 1); 
    x0 = rand(i, 1); 
    
    % Solución por gradiente conjugado 
    tic; 
    [NORMA1, iters1] = GradienteConjugadoVarianteB(A1, b, x0, 1e-6); 
    TIEMPO1var(i) = toc;

    % Solución por gradiente conjugado 
    tic 
    [NORMA2, iters2] = GradienteConjugadoVarianteB(A2, b, x0, 1e-6); 
    TIEMPO2var(i) = toc;

    % Solución por gradiente conjugado 
    tic; 
    [NORMA3, iters3] = GradienteConjugadoVarianteB(A3, b, x0, 1e-6); 
    TIEMPO3var(i) = toc;

end 

figure; 
plot(TIEMPO1); 
hold on; 
plot(TIEMPO2); 
plot(TIEMPO3); 
legend('Uniform', 'Clustered1', 'Clustered2'); 
title('TIEMPO vs n : Gradiente Conjugado VARIANTE'); 

% Vamos a comparar los tiempos para cada matriz y para cada método 
figure; 

subplot(3, 1, 1); 
plot(TIEMPO1); 
hold on; 
plot(TIEMPO1var); 
title('Comparando tiempo vs n para cada método y para cada matriz')

subplot(3, 1, 2); 
plot(TIEMPO2); 
hold on; 
plot(TIEMPO2var); 

subplot(3, 1, 3); 
plot(TIEMPO3); 
hold on; 
plot(TIEMPO3var); 


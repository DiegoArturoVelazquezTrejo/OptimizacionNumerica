d = 2; 
n = 2; 

x = rand(d, n); 

f = @(x) 10*((x(2)-x(1)^2)^2) + 1 - x(1)^2;
df = @(x)[ 40*x(1)^3-40*x(2)*x(1)-2*x(1); 20*(x(2)-x(1)^2)];
d2f = @(x) [120*x(1)^2-40*x(2)-2, -40*x(1); -40*x(1), 20];

N = 50; 
t = linspace(-5, 5, N); 
[X, Y] = meshgrid(t, t); 

fplot = @(x, y) f([x;y]); 
dfplot = @(x, y, j) df([x;y]); 
Z = arrayfun(fplot, X, Y);
G = zeros(N, N, d); 

for i= 1:N
    for j = 1:N 
        G(i,j,:) = dfplot(X(i, j), Y(i, j));
    end
end 
contour(X, Y, Z, 100); 
axis equal 
hold on; 
quiver(X, Y, G(:,:,1), G(:,:,2)); 
plot(x(1, :), x(2, :), 'r*');


MAXITER = 100; 
EPS = 1e-4; 
x0=[-5;-5];

%[xk, iteraciones, XSKS, FXKS] = BusquedaLinea(x0, f, df, df, MAXITER, EPS, "SR1"); 
delta_hat = 0.5; 
etha = 1e-2; % \in (0, 1/4)
[xk, iteraciones, XSKS, FXKS] = RegConf(f, df, x0, "dogleg", "BFGS", delta_hat, etha, EPS, MAXITER, d2f); 

plot(XSKS(1,:),XSKS(2,:),'k*--')
xk 
iteraciones
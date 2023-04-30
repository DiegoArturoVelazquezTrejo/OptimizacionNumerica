d = 2; 
n = 5; 

x = rand(d, n); 

f = @(z) sum(vecnorm(x-repmat(z, 1, n))); 
df = @(z) sum((repmat(z, 1, n) - x)./repmat(vecnorm(x-repmat(z, 1, n)), d, 1),2); 
N = 50; 
t = linspace(0, 1, N); 
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
contour(X, Y, Z, 20); 
axis equal 
hold on; 
quiver(X, Y, G(:,:,1), G(:,:,2)); 
plot(x(1, :), x(2, :), 'r*');


MAXITER = 1000; 
EPS = 1e-8; 
x0 = [1;1]; 

%[xk, iteraciones, XSKS, FXKS] = BusquedaLinea(x0, f, df, df, MAXITER, EPS, "SR1"); 
delta_hat = 0.25; 
etha = 1e-2; % \in (0, 1/4)
[xk, iteraciones, XSKS, FXKS] = RegConf(f, df, x0, "dogleg", "BFGS", delta_hat, etha, EPS, MAXITER, []); 

plot(XSKS(1,:),XSKS(2,:),'k*--')
xk 
iteraciones
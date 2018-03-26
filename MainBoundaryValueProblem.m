function MainBoundaryValueProblem
disp("Решение краевой задачи y'' + p(x)y' + q(x)y = f(x)");
disp("alpha_0 * y(a) + alpha_1 * y'(a) = A, |alpha_0| + |alpha_1| > 0");
disp("beta_0 * y(b) + beta_1 * y'(b) = B, |beta_0| + |beta_1| > 0");
disp("С помощью метода прогонки по основной сетке");
n = input('Число разбиений: ');
syms x y p q f ;
alpha = 0.4;
alpha0 = 1;
alpha1 = -2;
A = 1/alpha;
beta0 = 1;
beta1 = 0;
B = 1/(1 + alpha);
y(x) = 1 / (x^2 + alpha);
p(x) = - (x^2 + alpha);
q(x) = -2*x;
f(x) = 2 * (3*x^2 - alpha)/(x^2 + alpha)^3;
h = 1/n;
X = h.*(0:n)';
Y = zeros(n+1, 1);
Z = y(X);
P = p(X);
Q = q(X);
F = f(X);
a = 1 - h*P/2;
b = 1 + h*P/2;
c = 2 - (h^2) * Q;
p1 = alpha1 / (alpha1 - h*alpha0);
d1 = - A*h / (alpha1 - h*alpha0);
p2 = beta1 / (h*beta0 + beta1);
d2 = B*h / (h*beta0 + beta1);
m = zeros(n, 1);
m(1,1) = p1;
k = zeros(n, 1);
k(1,1) = d1;
for i = 2:n
    m(i, 1) = b(i) / (c(i) - a(i) * m(i-1,1));
    k(i, 1) = (a(i) * k(i-1,1) - (h^2) * F(i)) / (c(i) - a(i) * m(i-1, 1));
end
Y(n+1, 1) = (p2*k(n,1) + d2) / (1 - p2*m(n, 1));
for i = n:-1:1
    Y(i,1) = m(i,1)*Y(i+1, 1) + k(i, 1);
end
disp('Фактическая абсолютная погрешность:');
disp(vpa(abs(Z-Y), 6));
disp('Максимальное значение погрешности:');
disp(vpa(max(abs(Z-Y)), 6));
end


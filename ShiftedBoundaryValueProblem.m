function ShiftedBoundaryValueProblem
disp("Решение краевой задачи y'' + p(x)y' + q(x)y = f(x)");
disp("alpha_0 * y(a) + alpha_1 * y'(a) = A, |alpha_0| + |alpha_1| > 0");
disp("beta_0 * y(b) + beta_1 * y'(b) = B, |beta_0| + |beta_1| > 0");
disp("С помощью метода прогонки по сдвинутой сетке");
n = input('Число разбиений: ');
syms x y p q f ;
alpha = 0.4;
alpha_0 = 1;
alpha_1 = -2;
A = 1/alpha;
beta_0 = 1;
beta_1 = 0;
B = 1/(1 + alpha); 
y(x) = 1 / (x^2 + alpha);
p(x) = - (x^2 + alpha);
q(x) = -2*x;
f(x) = 2 * (3*x^2 - alpha)/((x^2 + alpha)^3);
h = 1/n;
X = (-h/2 + h.*(0:n+1))';
Y = zeros(n+2, 1);
Z = y(X);
P = p(X);
Q = q(X);
F = f(X);
a = ones(n+2, 1) - h*P/2;
b = ones(n+2, 1) + h*P/2;
c = 2*ones(n+2, 1) - h^2 * Q;
p_1 = (alpha_0*h + 2*alpha_1)/(2*alpha_1 - alpha_0*h);
d_1 = 2*A*h / (alpha_0*h - 2*alpha_1);
p_2 = (2*beta_1 - beta_0*h) / (beta_0*h + 2*beta_1);
d_2 = 2*B*h / (beta_0*h + 2*beta_1);
m = zeros(n+1, 1);
m(1,1) = p_1;
k = zeros(n+1, 1);
k(1,1) = d_1;
for i = 2:n+1
    m(i, 1) = b(i,1) / (c(i,1) - a(i,1)*m(i-1,1));
    k(i, 1) = (a(i,1)*k(i-1,1) - h^2 * F(i,1))/(c(i,1) - a(i,1)*m(i-1, 1));
end
Y(n+2, 1) = (p_2*k(n+1,1) + d_2) / (1 - p_2*m(n+1,1));
for i = n+1:-1:1
    Y(i,1) = m(i,1)*Y(i+1, 1) + k(i, 1);
end
disp('Фактическая абсолютная погрешность');
disp(vpa(abs(Z-Y), 6));
disp('Максимальное значение погрешности');
disp(vpa(max(abs(Z-Y)), 6));
end


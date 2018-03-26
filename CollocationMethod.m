function CollocationMethod
    disp('Решение краевой задачи');
    disp("-(4+x)/(5+2x) u'' + (x/2 - 1)u' + (1 + exp(x/2))u = 2 + x");
    disp("u'(-1) = u(1) + 2u'(1) = 0");
    disp('методом коллокаций');
    n = input('Число координатных функций: ');
    syms x u p q r f W P;
    p(x) = - (4 + x)/(5 + 2*x);
    q(x) = (x/2 - 1);
    r(x) = 1 + exp(x/2);
    f(x) = 2 + x;
    W(1,1) = x^2 + 2*x - 11;
    W(2, 1) = x^3 - 3*x + 2;
    for i = 3:n
        W(i,1) = (1-x^2)^2 * jacobiP(i-3,2,2,x);
    end
    A = zeros(n);
    F = zeros(n,1);
    T = cos((2* (1:n) -1)*pi / (2*n));
    for j = 1:n
        phi = p*diff(W(j,1), 2) + q*diff(W(j,1)) + r*W(j,1);
        A(:,j) = phi(T);
    end
    F(:,1) = f(T);
    c = A\F;
    disp('Расширенная матрица системы:');
    disp([A F]);
    disp('Число обусловленности матрицы A');
    disp(cond(A));
    disp('Коэффициенты разложения приближенного решения по координатным функциям');
    c = A\F;
    disp(c');
    disp('Значения приближенного решения в точках -0.5, 0, 0.5:');
    disp(vpa(subs(c'*W, x, -0.5 : 0.5 : 0.5), 8));
    u_n = c'*W;
    X = linspace(-1,1,100);
    Y = subs(u_n, x, X);
    plot(X,Y);
end


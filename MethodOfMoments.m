function MethodOfMoments
    disp('Решение краевой задачи');
    disp("-(4+x)/(5+2x) u'' + (x/2 - 1)u' + (1 + exp(x/2))u = 2 + x");
    disp("u'(-1) = u(1) + 2u'(1) = 0");
    disp('методом моментов');
    n = input('Число координатных функций: ');
    syms x u p q r f W P;
    p(x) = - (4 + x)/(5 + 2*x);
    q(x) = (x/2 - 1);
    r(x) = 1 + exp(x/2);
    f(x) = 2 + x;
    P(1,1) = jacobiP(0,0,0,x);
    P(2,1) = jacobiP(1,0,0,x);
    for i = 3:n
        P(i,1) = jacobiP(i-1,0,0,x);
    end
    W(1,1) = x^2 + 2*x - 11;
    W(2, 1) = x^3 - 3*x + 2;
    for i = 3:n
        W(i,1) = (1-x^2)^2 *jacobiP(i-3,2,2,x);
    end
    A = zeros(n);
    F = zeros(n,1);
    L = x.^(0 : 1);
    for k = 3 : 11
        L(1, k) = (2*(k - 1) - 1)/(k - 1) * x * L(1, k - 1) - ((k - 1) - 1)/(k - 1) * L(1, k - 2);
    end
    Leg = L(1, 11);
    T = solve(Leg);
    C = 2./((1 - T.*T).*(subs(diff(Leg, x), x, T)).^2);
    for i = 1:n
        for j = 1:n
            phi = (p*diff(W(j,1), 2) + q*diff(W(j,1)) + r*W(j,1))*P(i,1);
            A(i,j) = C'*phi(T);
            phi = f*P(i,1);
            F(i,1) = C'*phi(T);
        end
    end
    disp('Расширенная матрица системы:');
    disp([A F]);
    disp('Число обусловленности матрицы A');
    disp(cond(A));
    disp('Коэффициенты разложения приближенного решения по координатным функциям');
    c = A\F;
    disp(c');
    disp('Значения приближенного решения в точках -0.5, 0, 0.5:');
    u_n = c'*W;
    disp(vpa(subs(u_n, x, -0.5 : 0.5 : 0.5), 8));
    %u = dsolve('-(4 + x)/(5 + 2*x)*D2u + (x/2 - 1)*Du + (1 + exp(x/2))*u = 2 + x', 'Du(-1)=0, u(1)+2*Du(1)=0', 'x')
    %disp('Невязка:');
    %disp(vpa(subs(c'*W, x, -0.5 : 0.5 : 0.5) - subs(u, x, -0.5 : 0.5 : 0.5)));
    X = linspace(-1,1,100);
    Y = subs(u_n, x, X);
    plot(X,Y);
end


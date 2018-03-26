function SturmLiouville
    disp('Проблема собственных значений в задаче Штурма-Лиувилля');
    disp("-((kx+l)u')' + (k^2(1/(kx+l) - kx))u = ?u, k = 1.57894, l = 8.59453");
    disp('u(-1) = u(1) = 0');
    syms x p q y_n;
    global T C W G
    p(x) = 1.57894*x + 8.59453;
    q(x) = (1.57894)^2 *(1/(1.57894*x + 8.59453) - 1.57894*x);
    p_min =  p(-1);
    p_max = p(1);
    q_min =  q(1);
    q_max = q(-1);
    disp('Верхние и нижние оценки для первых двух собственных чисел и соответствующие им собственные функции оператора Штурма-Лиувилля:');
    disp(vpa((pi/2)^2 * p_min + q_min, 8));
    disp(vpa(subs(-diff(p_min*diff(cos(pi/2 *x))) + q_min*cos(pi/2 *x) - ....
        ((pi/2)^2 * p_min + q_min)*cos(pi/2 *x), 0), 13));
    disp(vpa((pi/2)^2 * p_max + q_max, 6));
    disp(vpa(subs(-diff(p_max*diff(cos(pi/2 *x))) + q_max*cos(pi/2 *x) - ....
        ((pi/2)^2 * p_max + q_max)*cos(pi/2 *x), 0), 13));
    disp('cos(pi/2 * x)');
    disp(vpa(pi^2* p_min + q_min, 6));
    disp(vpa(subs(-diff(p_min*diff(sin(pi*x))) + q_min*sin(pi*x) - ....
        ((pi)^2 * p_min + q_min)*sin(pi*x), 0), 13));
    disp(vpa(pi^2 * p_max + q_max, 6));
    disp(vpa(subs(-diff(p_max*diff(sin(pi*x))) + q_max*sin(pi*x) - ....
        ((pi)^2 * p_max + q_max)*sin(pi*x), 0), 13));
    disp('sin(pi*x)');
    disp('Приближения для первых двух собственных чисел по формуле:');
    phi = p*(diff(cos(pi/2 *x)))^2 + q*(cos(pi/2 *x))^2;
    disp(vpa(C'*phi(T), 8));
    phi = p*(diff(sin(pi*x)))^2 + q*(sin(pi*x))^2;
    disp(vpa(C'*phi(T), 8));
    e = sort(eig(G));
    disp(vpa(e(1)));
    disp(vpa(e(2)));
    epsilon = input('Точность: ');
    disp('Метод скалярных произведений');
        y = ones(size(G, 1), 1);
        z = G\y;
        lambda0 = 0;
        lambda1 = (z'*y)/(y'*y);
        while abs(lambda1 - lambda0) >= epsilon
            y = z;
            [~,I] = max(abs(y));
            y = y / y(I,1);
            z = G\y;
            lambda0 = lambda1;
            lambda1 = (z'*y)/(y'*y);
        end
        lambda1 = 1/lambda1;
        disp(vpa(lambda1));
end


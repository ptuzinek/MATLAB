format long
x = pi/2;
alfan = x-0.01*x:.0001:x+0.01*x;  
for i=1:length(alfan)
    [dm3,cm3] = matrix(alfan(i),3); %wartosc wyznacznika i wskaznika uwarunkowania
    [dm10,cm10] = matrix(alfan(i),10);
    [dm20,cm20] = matrix(alfan(i),20);
    d1(i) = dm3;
    d2(i) = dm10;
    d3(i) = dm20;
    c1(i) = cm3;
    c2(i) = cm10;
    c3(i) = cm20;   
end
% Zadanie 2
semilogy(alfan,d1,'b-'); %nieb. ci¹g³a bez markera
hold on;
title('Wykres 1. Wykres wyznaczników macierzy')
xlabel('argument macierzy alfa')
ylabel('wyznacznik macierzy A')
semilogy(alfan,d2,'r-'); %czerwona
hold on;
semilogy(alfan,d3,'k-'); %czarna
hold on;
grid on;
figure
semilogy(alfan,c1,'b-'); %niebieskie
hold on;
title('Wykres 2. Wykres wskaŸników uwarunkowania')
xlabel('argument macierzy alfa')
ylabel('wskaŸnik uwarunkowania macierzy A')
semilogy(alfan,c2,'r-'); %czerwona
hold on;
semilogy(alfan,c3,'k-'); %czarna
hold on;
grid on;
% Zadanie 4
for  N = [3 10 20]
x1 = [1:N]';
t = pi/2;
alfan = t-0.01*t:.0001:t+0.01*t;
A = macierz(alfan, N);
b = macierz(alfan, N) * x1;
CB = CB_rozw2(A,b)
end
% Zadanie 5
p=2;
for N = [3 10 20]
    k = 0;
    t = pi/2;
    xalfan = t-0.01*t:.0001:t+0.01*t;
    x=1:N;
    w_kwadrat = t-0.01*t:.0001:t+0.01*t;
    w_max = t-0.01*t:.0001:t+0.01*t;
    w_kwadratpseudo = t-0.01*t:.0001:t+0.01*t;
    w_maxpseudo = t-0.01*t:.0001:t+0.01*t;
    for  alfan = t-0.01*t:.0001:t+0.01*t;
        k=k+1;
        A = macierz(alfan, N);
        b = A*x';
        X1 = CB_rozw2(A,b);
        detX1 = X1-x';
        X2 = A\b;
        detX2 = X2-x';
        kwadrat = norm(detX1,2)/norm(x,2);
        max = norm(detX1,inf)/norm(x,inf);
        kwadratpseudo = norm(detX2,2)/norm(x,2);
        maxpseudo = norm(detX2,inf)/norm(x,inf);
        w_kwadrat(1,k) = kwadrat;
        w_max(1,k) = max;
        w_kwadratpseudo(1,k)= kwadratpseudo;
        w_maxpseudo(1,k)= maxpseudo;  
    end
p = p+1;
figure
semilogy(xalfan, w_kwadrat,'b-', xalfan, w_kwadratpseudo,'g-');
str = ['Wykres ',num2str(p),'. Wykres b³êdów kwadratowych rozwi¹zañ równania dla N=',num2str(N)];
title(str);
grid on;
p=p+1;
figure
semilogy(xalfan, w_max,'b', xalfan, w_maxpseudo,'g');
str = ['Wykres ',num2str(p),'. Wykres b³êdów maksymalnych rozwi¹zañ równania dla N=',num2str(N)];
title(str);
grid on;
end
% Zadanie 3
function  x  = CB_rozw2( A, b )
%   Funkcja wyznacza rozwiazania ukladu rownan A*x = b poprzez rozklad CB
[N,m]=size(A);
L=zeros(N);
% tworzenie sumy elementów macierzy L
for n = 1:N                 
    l1=0;
    for i = 1:(N-1)
        l1 = l1 + ((L(n,i))^2); 
    end
    L(n,n) = sqrt(A(n,n) - l1);
    for v = n+1:N
        l2=0;
        for i = 1:(n-1) 
            l2 = l2 + L(v,i) * L(n,i);   
        end
        L(v,n) = (A(v,n)- l2)/L(n,n);
    end
end
y = inv(L)*b;
x = inv(L')*y;
end
% Utworzenie macierzy z zadania 1
function  A  = macierz( alfan, N )
format long
A = ones(N);
A(1,1) = (f(alfan))^2;
i=2;
x=0;
    for n = 2:N
    A(1,n) = ((-1)^(n+1))*2*f(alfan);
    A(n,1) = ((-1)^(n+1))*2*f(alfan);
    x=x+1;
    for i=n:N
        if  rem(x,2) == 0
             A(n,i) = ((-1)^(i+1))*4*n;
             A(i,n) = ((-1)^(i+1))*4*n;
        else 
            A(n,i) = ((-1)^(i))*4*n;
            A(i,n) = ((-1)^(i))*4*n;
        end
    end
    end

end
% Zadanie 2
function [y,w] = matrix( alfan, N )
format long
A = ones(N);
A(1,1) = (f(alfan))^2;
i=2;
x=0;
    for n = 2:N
    A(1,n) = ((-1)^(n+1))*2*f(alfan);
    A(n,1) = ((-1)^(n+1))*2*f(alfan);
    x=x+1;
    for i=n:N
        if  rem(x,2) == 0
             A(n,i) = ((-1)^(i+1))*4*n;
             A(i,n) = ((-1)^(i+1))*4*n;
        else 
            A(n,i) = ((-1)^(i))*4*n;
            A(i,n) = ((-1)^(i))*4*n;
        end
    end
    end
   
y= det(A);
w= cond(A);
end
% Zapewnienie zera
function  y  = f( a )
format long
y = cos(a)/a;
if a == pi/2
    y = 0;
end
end

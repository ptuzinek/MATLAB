% zadanie 1
clear all
s=1;
for N = [5 10 15]
    K=N;
    n=1:N;
    x = -1 + 2*(n-1)/(N-1)
    for i = 1:N
        f(i) = cos(pi*x(i))*exp(x(i))
    end
    figure
    plot(x,f,'ro')
    str=['Wykres ',num2str(s),'. Wykres funkcji dla N = ',num2str(N)]
    title(str)
    xlabel('xn')
    ylabel('f(xn)')
    zlabel('B��d �redniokwadratowy')
    %legend('funkcja f')
    s=s+1;
end
% zadanie 2
figure
zadanie22(22,20,1,s)
s=s+1;
figure
zadanie22(22,20,2,s)
s=s+1;
figure
zadanie22(10,5,1,s)
s=s+1;
figure
zadanie22(10,5,2,s)
% zadanie 3
clear all
m = input(' Wybierz metod� rozwi�zania uk�adu r�wna�. 1 - rozk�ad CB, 2 - pseudo-odwrotno�� macierzy ')
xg= linspace(-1,1,1000);
r = 5:50;
s=8;
for N = 5:50
    for K =5:N
        ap = aproksymacja(N,K,m);
        for i = 1:1000
            y(i) = cos(pi*xg(i))*exp(xg(i));
        end
    dety = ap-y';
    sk(N,K) = norm(dety,2)/norm(y,2); % zapis b�ed�w kolejnych punkt�w dla danego N w jednym wierszu
    max(N,K) = norm(dety,inf)/norm(y,inf);
    end
end
% WYKRES B��DU �REDNIOKWADRATOWEGO
figure
surf(log10(sk))
str=['Wykres ',num2str(s),'. Wykres b��d�w �redniokwadratowych rozwi�za� r�wnania '];
title(str)
xlabel('N')
ylabel('K')
zlabel('B��d �redniokwadratowy')
s=s+1;
% WYKRES B��DU MAKSYMALNEGO
figure
surf(log10(max))
str=['Wykres ',num2str(s),'. Wykres b��d�w maksymalnych rozwi�za� r�wnania '];
title(str)
xlabel('N')
ylabel('K')
zlabel('B��d maksymalny')

% zadanie 4       

m = input(' Wybierz metod� rozwi�zania uk�adu r�wna�. 1 - rozk�ad CB, 2 - pseudo-odwrotno�� macierzy ');
r = 5:50;
w=10;
s=0;
wynik1=zeros(50,50,6);
wynik2=zeros(50,50,6);
xg= linspace(-1,1,1000);
for i = 1:1000
            f(i) = cos(pi*xg(i))*exp(xg(i));
end
for si = .1:0.3:1.6
s=s+1;
for N = 5:50
    dy= si*randn(1,N);
    for K =5:N
        apz = zaproksymacja(N,K,m,si);
     dety = apz-f';
     sk(N,K) = norm(dety,2)/norm(f,2); % zapis b�ed�w kolejnych punkt�w dla danego N w jednym wierszu
     max(N,K) = norm(dety,inf)/norm(f,inf);
    end
end
wynik1(:,:,s) = sk;
wynik2(:,:,s) = max;
end
sigma = .1:0.3:1.6;
% WYKRES ZALE�NO�CI B��D�W OD SIGMY
for i = 1:6
    si = sigma(i);
    figure
    surf(log10(wynik1(:,:,i)))
    str=['Wykres ',num2str(w),'. Wykres b��d�w �redniokwadratowych rozwi�za� r�wnania dla zaburzonych danych, si = ',num2str(si)];
    title(str)
    xlabel('N')
    ylabel('K')
    zlabel('B��d �redniokwadratowy')
    w=w+1;
    figure
    surf(log10(wynik2(:,:,i)))
    str=['Wykres ',num2str(w),'. Wykres b��d�w maksymalnych rozwi�za� r�wnania dla zaburzonych danych, si = ',num2str(si)];
    title(str)
    xlabel('N')
    ylabel('K')
    zlabel('B��d maksymalny') 
    w=w+1;
end

function [ apz ] = zaproksymacja(N,K,m,si )

n=1:N;
xn = -1 + 2*(n-1)/(N-1);
k=1:K;
xk = -1 + 2*(k-1)/(K-1);
for i = 1:N   
    for j = 1:K
        G(i,j) = (4/sqrt(2*pi))*exp(-16*(xn(i)-xk(j))^2);
    end
    y(i) = cos(pi*xn(i))*exp(xn(i));
end
    dy = si*randn(1,N);
       yz = y + dy;

       A=G'*G;
       bz=G'*yz';
       pz=A\bz;
if m == 1
    pz=CB_rozw2(A,bz);
elseif m == 2
    pz=A\bz;
else
    a=('B��d, wybierz metod� rozwi�zania uk�adu')
end

xg= linspace(-1,1,1000); 
        for i = 1:1000
            for j = 1:K
                G1(i,j) = (4/sqrt(2*pi))*exp(-16*(xg(i)-xk(j))^2);
            end
        end
apz=G1*pz;
end





function [ ap ] = aproksymacja( N, K, m )
%   funcja przybli�ana - y, macierz przekszta�caj�ca - G
%   Dostje N i K, wylicza wektory xn, xk, macierz Fi oraz wektor p
%   Nast�pnie wylicza warto�� funkcji aproksymuj�cej dla wielu punkt�w
%   Rysuje wykres je�li odkomentuje si� polecenia
%   Wybierz metod� rozwi�zania uk�adu r�wna�. 1 - rozk�ad CB, 2 - pseudo-odwrotno�� macierzy "\"
n=1:N;
xn = -1 + 2*(n-1)/(N-1);
k=1:K;
xk = -1 + 2*(k-1)/(K-1);
for i = 1:N   
    for j = 1:K
        G(i,j) = (4/sqrt(2*pi))*exp(-16*(xn(i)-xk(j))^2);
    end
    y(i) = cos(pi*xn(i))*exp(xn(i));
end
A=G'*G;
b=G'*y';
if m == 1
    p=CB_rozw2(A,b);
elseif m == 2
    p=A\b;
else
    a=('B��d, wybierz metod� rozwi�zania uk�adu')
end
xg= linspace(-1,1,1000); 
        for i = 1:1000
            for j = 1:K
                G1(i,j) = (4/sqrt(2*pi))*exp(-16*(xg(i)-xk(j))^2);
            end
        end
ap=G1*p;
end


function [ ap ] = zadanie22( N, K, m,s )
%   funcja przybli�ana - y, macierz przekszta�caj�ca - G
%   Dostje N i K, wylicza wektory xn, xk, macierz Fi oraz wektor p
%   Nast�pnie wylicza warto�� funkcji aproksymuj�cej dla wielu punkt�w
%   Rysuje wykres je�li odkomentuje si� polecenia
%   Wybierz metod� rozwi�zania uk�adu r�wna�. 1 - rozk�ad CB, 2 - pseudo-odwrotno�� macierzy "\"
n=1:N;
xn = -1 + 2*(n-1)/(N-1);
k=1:K;
xk = -1 + 2*(k-1)/(K-1);
for i = 1:N   
    for j = 1:K
        G(i,j) = (4/sqrt(2*pi))*exp(-16*(xn(i)-xk(j))^2);
    end
    y(i) = cos(pi*xn(i))*exp(xn(i));
end
A=G'*G;
b=G'*y';
if m == 1
    p=CB_rozw2(A,b);
elseif m == 2
    p=A\b;
else
    a=('B��d, wybierz metod� rozwi�zania uk�adu')
end
xg= linspace(-1,1,1000); 
        for i = 1:1000
            for j = 1:K
                G1(i,j) = (4/sqrt(2*pi))*exp(-16*(xg(i)-xk(j))^2);
            end
        end
ap=G1*p;
plot(xg,ap, xn,y,'ro')
str=['Wykres ',num2str(s),'. Wykres w�z��w funkcji f oraz aproksymacji'];
title(str)
xlabel('xn')
ylabel('Warto�� funkcji')
legend('Aproksymacja','W�z�y')
end

function  x  = CB_rozw2( A, b )
%   Funkcja wyznacza rozwiazania ukladu rownan A*x = b poprzez rozklad CB
[N,m]=size(A);
L=zeros(N);
% tworzenie sumy element�w macierzy L
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
% x=x';
end
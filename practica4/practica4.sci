funcprot(0);
exec('../lineales.sci',-1);

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 1.1 (4.2)==")
disp('(a) Como  U es inversible, su determinante no es cero. Al ser triangular superior, el producto de los elementos de la diagonal es igual al determinante. Sigue que ningún elemento en la diagonal es igual a cero. Teniendo esto en cuenta y planteando al ecuacion UX=I, se puede probar que mirando la columna i-ésima de I, los ultimos (i-1) elementos de X tienen que ser 0. Lo que significa que X es triangular superior, más aún X es la inversa de U y por únicidad de la inversa el ejercicio queda demostrado.');
disp('(b) ');
disp('(c) Se puede probar facilmente aplicando la formula de multiplicación de matrices.');

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 1.5 (4.2)==");
disp("Por teorema sabemos que existe L y U tal que A=LU => det(A)=det(L)det(U). det(A)!=0 por ser A no singular. Por lo que det(U)!=0. Como U es triangular, el producto de los elementos de la diagonal es igual al determinante, ningún elemento puede ser cero.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 1.8 (4.2)==");
disp("Podemos ver que existe descomposición LU. Por lo que es posible resolver el sistema Ax=e_i para todo i. Resolviendo estos n sistemas y concatenando los resultados, tenemos la inversa.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 1.23 (4.2)==");
disp("La division por u_jj se debe realizar al final.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 1.36 (4.2)==");
disp("La matriz es real, simétrica y definida positiva. Queremos x/ Ax=b => (LL^T)x=b => Lz=b, L^Tx=z.");
A=[0.05 0.07 0.06 0.05; 0.07 0.10 0.08 0.07; 0.06 0.08 0.10 0.09; 0.05 0.07 0.09 0.10];
b=[0.23; 0.32; 0.33; 0.31];
L = Cholesky(A);
z = SustitucionProgresiva(L, b);
disp(z, "z = ");
x = SustitucionRegresiva(L', z);
disp(x, "x = ");

disp(norm(A*x-b), "El error de la solución es ");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 2==");
disp("Usamos la función Doolittle.");
A=[1 2  3  4;
   1 4  9  16;
   1 8  27 64;
   1 16 81 256];
b=[2; 10; 44; 190]
[L, U] = Doolittle(A);
z = SustitucionProgresiva(L, b);
x = SustitucionRegresiva(U, z);
disp(x, "x = ");
disp(norm(A*x-b), "El error de la solución es ");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 3==");
A=[16 -12 8 -16;
   -12 18 6 9;
   8 -6 5 -10;
   -16 9 -10 46];
B=[4 12 -16;
   12 37 -43;
   -16 -43 98];
//disp(chol(A)-Cholesky(A));
disp("La matriz A no es simetrica (ni definida positiva) por lo que chol() falla. Sin embargo nuestra función nos da un resultado aproximado:");
L = Cholesky(A);
disp(L*L', "L*L'' = ");
disp(norm(chol(B)'-Cholesky(B)), "Para la matriz B, no hay diferencia. norm(chol(B)'', Cholesky(B)) = ");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 4==");
disp("Ya se hizo en el TP2.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 5==");
disp("Se elaboró la función Tridiagonal(A, b).");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 6.20 (4.3)==");
disp("TODO.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 6.22 (4.3)==");
disp("Parecido al TP2.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 6.53 (4.3)==");
function x=Solve(A, b)
  _, P, L, U = Gauss(A,b);
  z = SustitucionProgresiva(L, P*b);
  x = SustitucionRegresiva(U, z);
endfunction
disp("P^-1 = P^T, PA=LU => A^T = U^T*L^T*P. y^T*A=c^T => A^T*y=c => U^T*L^T*P*y=c, el cual es fácilmente resolvible.")
function y=Solve2(A, c)
  _, P, L, U = Gauss(A,b);
  z = SustitucionProgresiva(U', c);
  u = SustitucionRegresiva(L', z);
  y = P'*u;
endfunction

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 7.2 (4.4)==");
disp("La demostracion salen usando las definiciones. La igualdad puede darse cuando x=e_i.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 7.7 (4.4)==");
disp("Solo hay que probar las propiedades.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 7.11 (4.4)==");
disp("Hay que usar la definicion aplicar las definiciones e ir despejando. Luego proponer un u y demostrar que ||Au|| es la norma.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 7.15 (4.4)==");
disp("Trivial.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 7.21 (4.4)==");
disp("Sabemos que ||A||=sup{||Ax||/||x||, x!=0}. ||A||=10. Sea x un vector cualquiera. 810=||A||>=||Ax||/||x|| => ||Ax||<=10||x||<=10. Vemos que la igualdad se cumple para x=[1;1;0] (||Ax||=10). ");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 7.40 (4.4)==");
disp("Trivial aplicando definicion y propiedades de la norma.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 7.49 (4.4)==");
disp("Trivial aplicando definición y propiedades de la norma e inversión de matrices.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 8.2 (4.6)==");
disp("Tomando la norma || ||_inf, podemos demostrar con la definicion que ||I-A||_inf<1. Por Teorema 1, el método de richardson (Q=I) converge a la solución.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 8.4 (4.6)==");
disp("Es el metodo de Gauss-Seidel tomando en cuenta que a_ii=0. Por Teorema sabemos que este método converge para matrizes diagonalmente dominante.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 8.7 (4.6)==");
disp("?? TODO");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 8.22 (4.6)==");
disp("Trivial despejando al fórmula.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 8.24 (4.6)==");
disp("Primera parte sale haciendo r^(k+1)=b-Ax^k y reemplazando x^k según la ecuación (3) despejada.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 8.27 (4.6)==");
disp("");
A=[3 1 1; 1 3 -1; 3 1 -5];
b=[5; 3; 1];
x=GaussSeidel(A, b, 30);
disp(x, "x = ");
disp(norm(A*x-b), "error = ");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 8.35 (4.6)==");
disp("Probar que (1-d)||x^(k)-x||<=d||x^(k) - x^(k-1)||, 1ero aplicar distributiva, usar el ''truquito'' -||A|| ||B||<= -||AB|| y la desigualdad triangular de la norma.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 9.a==");
disp("--> 1er sistema:");
disp("Se reorganizan las filas.");
A=[1 -1 -1; 1 -1 2; 0 2 4];
b=[0.375; 0; 0];
I=eye(3, 3);
[D, C_L, C_U] = Decompose(A);
Q=D;
disp(norm(I-inv(Q)*A), "Con Jacobi tenemos ||I-Q^-1 A|| = ");
disp("No puedo asegurarlo pues la matriz no es diagonalmente dominante, ni ||I-Q^-1 A||<1.");

disp("--> 2do sistema:");
disp("Supongo que la primera ecuación es x_1 - x_2 = 0");
A=[1 -2 0; 0 -1 1.1; -1 2 -1];
b=[0; 1; 0];
I=eye(3, 3);
[D, C_L, C_U] = Decompose(A);
Q=D;
disp(norm(I-inv(Q)*A), "Con Jacobi tenemos ||I-Q^-1 A|| = ");
disp("No puedo asegurarlo pues la matriz no es diagonalmente dominante, ni ||I-Q^-1 A||<1.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 9.b==");
A=[1 -1 -1; 1 -1 2; 0 2 4];
b=[0.375; 0; 0];
I=eye(3, 3);
[D, C_L, C_U] = Decompose(A);
Q=-C_L;
disp("Con Gauss-Seidel tenemos Q=D-C_L no es invertible.");
disp("No puedo asegurarlo pues la matriz no es diagonalmente dominante, ni ||I-Q^-1 A||<1.");

disp("--> 2do sistema:");
disp("Supongo que la primera ecuación es x_1 - x_2 = 0");
A=[
  1 -1 0;
  -1 2 -1;
  0 -1 1.1;
];
b=[0; 1; 0];
I=eye(3, 3);
[D, C_L, C_U] = Decompose(A);
Q=D-C_L;
disp(norm(I-inv(Q)*A), "Con Gauss-Seidel tenemos ||I-Q^-1 A|| = ");
disp("No puedo asegurarlo pues la matriz no es diagonalmente dominante, ni ||I-Q^-1 A||<1.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 9.c==");
disp("--> 1er sistema:");
disp("Se reorganizan las filas.");
A=[1 -1 -1; 1 -1 2; 0 2 4];
b=[0.375; 0; 0];
x=Jacobi(A, b, 30);
disp(norm(A*x-b), "error con Jacobi = ");
x=GaussSeidel(A, b, 30);
disp(norm(A*x-b), "error con GaussSeidel = ");

disp("--> 2do sistema:");
A=[
  1 -1 0;
  -1 2 -1;
  0 -1 1.1;
];
b=[0; 1; 0];
// disp(inv(A)*b); //Para espiar la solucion jeje
x=Jacobi(A, b, 100);
disp(norm(A*x-b), "error con Jacobi = ");
x=GaussSeidel(A, b, 100);
disp(norm(A*x-b), "error con GaussSeidel = ");
disp("El método no converge! (preguntar)");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 10==");
disp("Por el método de Gauss-Seidel tenemos que G=(D-C_L)^-1 C_U. Como ejemplo veamos cuando n=6.");
n=6;
A = zeros(n, n);
for i=1:n
  A(i, i) = 2;
  if i>1 then A(i, i-1) = -1; end
  if i<n then A(i, i+1) = -1; end
end
disp(A, "A = ");
[D, C_L, C_U] = Decompose(A);
disp(inv(D-C_L), "(D-C_L)^-1 = ");
disp("Se demuestra que ((D-C_L)^-1)_ij = 2^(i-j-1) si i>=j.");
for i=1:n
  for j=2:i+1
    X(i, j) = 2^(j-i-2);
  end
end
//disp(X);
disp(inv(D-C_L)*C_U, "G = ");
disp("Luego se demuestra que G_(i+1, j) = ((D-C_L)^-1)_ij, lo cual da una forma general para G.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 11==");
A=[10 1 2 3 4;
    1 9 -1 2 -3;
    2 -1 7 3 -5;
    3 2 3 12 -1;
    4 -3 -5 -1 15];
b = [12; -27; 14; -17; 12];
x=Jacobi(A, b, 50);
disp(norm(A*x-b), "error con Jacobi = ");
x=GaussSeidel(A, b, 50);
disp(norm(A*x-b), "error con GaussSeidel = ");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 12.1 (4.7)==");
disp("Trabajar y trabajar la expresión (gradiente(q))_i hasta llegar a (2(Ax-b))_i.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 12.3 (4.7)==");
disp("Se prueba demostrando que el producto escalar es 0. Para eso escribimos v_{k+1} como b-Ax_{k+1}, expandiendo el x_{k+1} también.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 13==");
disp("Se aproxima la solución del ejercicio 11.");
A=[10 1 2 3 4;
    1 9 -1 2 -3;
    2 -1 7 3 -5;
    3 2 3 12 -1;
    4 -3 -5 -1 15];
b = [12; -27; 14; -17; 12];
x=DescensoRapido(A, b, 325);
disp(norm(A*x-b), "error con DescensoRapido = ");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 14.a==");
disp("TODO");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 14.b==");
A=[1 0 0;
  -1 0 1;
  -1 -1 2];
[autovectores, autovalores] = spec(A);
disp(spec(A), "(a)");
A=[1 0 0;
  -0.1 0 0.1;
  -0.1 -0.1 2];
disp(spec(A), "(b)");
A=[4.75 2.25 -0.25;
  2.25 4.75 1.25;
  -0.25 1.25 4.75];
disp(spec(A), "(c)");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 15==");
disp("--> Ej 10");
disp("TODO");

disp("--> Ej 11");
disp("TODO Aplicar a la G del 11!");
A=[10 1 2 3 4;
    1 9 -1 2 -3;
    2 -1 7 3 -5;
    3 2 3 12 -1;
    4 -3 -5 -1 15];
b = [12; -27; 14; -17; 12];
[x, r] = MetodoPotencia(A, 500);
disp(x, "x =");
disp(r, "r =");
disp(norm(A*x-r*x), "error1 = ");
disp(abs(max(abs(spec(A)))-r), "error2 = ");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 16==");
disp("El método de la potencia inversa sirve para encontrar el autovalor r de A más pequeño. Utiliza el hecho de que 1/r es autovalor de A^-1 y es el mayor. Por lo que se podría sacar la inversa A^-1 y hacer el método de la potencia común sobre él, sin embargo esto trae muchos problemas de redondeo. Es mucho mejor utilizar Gauss para resolver Ax^(k+1)=x^(k). Y se puede aprovechar que el proceso de eliminación sólo debe ser realizado una vez.");
A=[1 0 0;
  -1 0 1;
  -1 -1 2];
A=[1 0 -0.5;
  -1 0 1;
  -1 -1 2];
  
disp(spec(A), "Autovalores = ");
[x, r] = MetodoPotenciaInversa(A, 500);
disp(r, "r =");
disp(x, "x =");
disp(norm(A*x-r*x), "error1 = ");
disp(abs(min(abs(spec(A)))-r), "error2 = ");




funcprot(0);
exec('../lineales.sci',-1);

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 1.1==")
disp('(a) Como  U es inversible, su determinante no es cero. Al ser triangular superior, el producto de los elementos de la diagonal es igual al determinante. Sigue que ningún elemento en la diagonal es igual a cero. Teniendo esto en cuenta y planteando al ecuacion UX=I, se puede probar que mirando la columna i-ésima de I, los ultimos (i-1) elementos de X tienen que ser 0. Lo que significa que X es triangular superior, más aún X es la inversa de U y por únicidad de la inversa el ejercicio queda demostrado.');
disp('(b) ');
disp('(c) Se puede probar facilmente aplicando la formula de multiplicación de matrices.');

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 1.5==");
disp("Por teorema sabemos que existe L y U tal que A=LU => det(A)=det(L)det(U). det(A)!=0 por ser A no singular. Por lo que det(U)!=0. Como U es triangular, el producto de los elementos de la diagonal es igual al determinante, ningún elemento puede ser cero.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 1.8==");
disp("Podemos ver que existe descomposición LU. Por lo que es posible resolver el sistema Ax=e_i para todo i. Resolviendo estos n sistemas y concatenando los resultados, tenemos la inversa.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 1.23==");
disp("La division por u_jj se debe realizar al final.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 1.36==");
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
disp("==Ejercicio 6.20==");
disp("TODO.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 6.22==");
disp("Parecido al TP2.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 6.53==");
function x=Solve(A, b)
  U, _, L, P = Gauss(A,b);
  z = SustitucionProgresiva(L, P*b);
  x = SustitucionRegresiva(U, z);
endfunction
disp("P^-1 = P^T, PA=LU => A^T = U^T*L^T*P. y^T*A=c^T => A^T*y=c => U^T*L^T*P*y=c, el cual es fácilmente resolvible.")
function y=Solve2(A, c)
  U, _, L, P = Gauss(A,b);
  z = SustitucionProgresiva(U', c);
  u = SustitucionRegresiva(L', z);
  y = P'*u;
endfunction

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 7.2==");
disp("La demostracion salen usando las definiciones. La igualdad puede darse cuando x=e_i.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 7.7==");
disp("Solo hay que probar las propiedades.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 7.11==");
disp("Hay que usar la definicion aplicar las definiciones e ir despejando. Luego proponer un u y demostrar que ||Au|| es la norma.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 7.15==");
disp("Trivial.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 7.21==");
disp("Sabemos que ||A||=sup{||Ax||/||x||, x!=0}. ||A||=8. Sea x un vector cualquiera. 8=||A||>=||Ax||/||x|| => ||Ax||<=8||x||<=8. Vemos que la igualdad se cumple para x=[1;1;0] (||Ax||=8). ");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 7.40==");
disp("Trivial aplicando definicion y propiedades de la norma.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 7.49==");
disp("Trivial aplicando definición y propiedades de la norma e inversión de matrices.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 8.2==");
disp("Tomando la norma || ||_inf, podemos demostrar con la definicion que ||I-A||_inf<1. Por Teorema 1, el método de richardson (Q=I) converge a la solución.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 8.4==");
disp("Es el metodo de Gauss-Seidel tomando en cuenta que a_ii=0. Por Teorema sabemos que este método converge para matrizes diagonalmente dominante.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 8.7==");
disp("??");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 8.22==");
disp("Trivial despejando al fórmula.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 8.24==");
disp("Primera parte sale haciendo r^(k+1)=b-Ax^k y reemplazando x^k según la ecuación (3) despejada.");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 8.27==");
disp("");
A=[3 1 1; 1 3 -1; 3 1 -5];
b=[5; 3; 1];
x=GaussSeidel(A, b, 30);
disp(x, "x = ");
disp(norm(A*x-b), "error = ");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 8.35==");
disp("TODO");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 9.a==");
disp("No puedo asegurarlo pues las matrices no son diagonalmente dominantes.");
disp("Supongo que la primera ecuación del segundo sistema es x_1 - x_2 = 0");
/////////////////////////////////////////////////////////////////
disp("==Ejercicio 9.b==");
disp("Lo mismo que 9.a?");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 9.c==");
disp("Se reorganizan las filas para hacer las matrices lo mas ''diagonalmente dominantes'' posibles.");
A=[1 -1 -1; 1 -1 2; 0 2 4];
b=[0.375; 0; 0];
x=Jacobi(A, b, 30);
disp(norm(A*x-b), "error con Jacobi = ");
x=GaussSeidel(A, b, 30);
disp(norm(A*x-b), "error con GaussSeidel = ");

A=[1 -2 0; 0 -1 1.1; -1 2 -1];
b=[0; 1; 0];
// disp(inv(A)*b); //Para espiar la solucion jeje
x=Jacobi(A, b, 50);
disp(norm(A*x-b), "error con Jacobi = ");
x=GaussSeidel(A, b, 50);
disp(norm(A*x-b), "error con GaussSeidel = ");
disp("TODO: arreglar, no convergen!");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 10==");
disp("TODO");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 11.==");
disp("TODO");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 12.1==");
disp("TODO");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 12.3==");
disp("TODO");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 13==");
disp("TODO");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 14.a==");
disp("TODO");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 14.b==");
disp("TODO");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 15==");
disp("TODO");

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 16==");
disp("TODO");





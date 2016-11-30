//Alumno: Martín Villagra
funcprot(0);

A=[3   2    0.5 ;
   2   3    0.75;
   0.5 0.75 3   ];
b=[6; 7.75; 4.75];
sol=[0.5; 2; 1];

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 2.b)==");

function [x, it, tasa]=Jacobi(A, b, tol, sol, x)
  n = length(b);
  if ~exists("x","local") then
    for i=1:n x(i)=0; end
  end
  tasa = -1;
  for k=1:10000
    if norm(x-sol)<tol then //Reviso si la iteración anterior fue suficiente
        it = k-1;
        return;
    end
    for i=1:n
      sum = 0;
      for j=1:n
        if j<>i then
          sum = sum + A(i, j)*x(j);
        end
      end
      u(i) = (b(i) - sum) / A(i, i);
    end
    tasa = norm(u-sol,%inf)/norm(x-sol, %inf);
    for i=1:n
      x(i) = u(i);
    end
  end
  it = -1; // No converge
endfunction

function [x, it, tasa]=GaussSeidel(A, b, tol, sol, x)
  n = length(b);
  if ~exists("x","local") then
    for i=1:n x(i)=0; end
  end
  tasa = -1;
  for k=1:10000
    if norm(x-sol)<tol then //Reviso si la iteración anterior fue suficiente
        it = k-1;
        return;
    end
    for i=1:n
      oldx(i) = x(i);
    end
    for i=1:n
      sum = 0;
      for j=1:n
        if j<>i then
          sum = sum + A(i, j)*x(j);
        end
      end
      x(i) = (b(i) - sum) / A(i, i);
    end
    tasa = norm(x-sol,%inf)/norm(oldx-sol, %inf);
  end
  it = -1; // No converge!
endfunction

function C=dot(A, B)
  C=A'*B;
endfunction

function [x, it, tasa]=DescensoRapido(A, b, tol, sol, x)
  n = length(b);
  if ~exists("x","local") then
    for i=1:n x(i)=0; end
  end
  tasa = -1;
  for k=1:10000
    if norm(x-sol)<tol then //Reviso si la iteración anterior fue suficiente
        it = k-1;
        return;
    end
    v = b - A*x;
    div = dot(v, A*v);
    if abs(div)<%eps then break; end
    t = dot(v, v) / div;
    nx = x + t*v;
    tasa = norm(nx-sol,%inf)/norm(x-sol, %inf);
    x = nx;
 end
  it = -1; // No converge!
endfunction

function its=IteracionesNecesarias(A, b, tol, sol)
    [_, jacob]=Jacobi(A, b, tol, sol);
    [_, seidel]=GaussSeidel(A, b, tol, sol);
    [_, descenso]=DescensoRapido(A, b, tol, sol);
    its=[jacob seidel descenso];
endfunction

disp(IteracionesNecesarias(A, b, 10^-7, sol), "El vector solicitado es = ");
// El vector solicitado es = 68. 22. 42. 

format(23);
/////////////////////////////////////////////////////////////////
disp("==Ejercicio 2.c)i)==");
function [D, C_L, C_U]=Decompose(A)
  n = size(A, 'r');
  C_L = 0*A;
  for i=1:n
    D(i, i) = A(i, i);
    for j=1:i-1
        C_L(i, j) = -A(i, j);
    end
  end
  C_U = -(A - D + C_L);
endfunction

function [r, x]=MetodoPotencia(A, it, sol, x)
  n = size(A, 'r');
  if ~exists("x","local") then
    for i=1:n x(i)=1; end
  end
  for k=1:it
    y = A*x;
    r = y(1)/x(1);
    x = y/norm(y);
 end
endfunction

function p=RadioEspectral(A)
    p=max(abs(spec(A)))
endfunction

[D, C_L, C_U]=Decompose(A);

disp("--> Jacobi");
G=inv(D)*(C_L+C_U);
disp(RadioEspectral(G), "Tasa de convergencia utilizando spec =");
// Tasa de convergencia utilizando spec = 0.77851351303733784537

disp(abs(MetodoPotencia(G, 50)), "Tasa de convergencia utilizando el método de la potencia =");
// Tasa de convergencia utilizando el método de la potencia = 0.77855980886291342724

disp("--> Gauss-Seidel");
G=inv(D-C_L)*C_U;
disp(RadioEspectral(G), "Tasa de convergencia utilizando spec =");
// Tasa de convergencia utilizando spec = 0.44444444444444441977

disp(abs(MetodoPotencia(G, 50)), "Tasa de convergencia utilizando el método de la potencia =");
// Tasa de convergencia utilizando el método de la potencia = 0.44444444444444441977

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 2.c)ii)==");
disp("Se adaptó la función de los métodos para que devuelvan la tasa de convergencia de la última iteración.");

disp("--> Jacobi");
[_, it, tasa]=Jacobi(A, b, 10^-7, sol);
disp(tasa, "Tasa de convergencia aproximada =");
// Tasa de convergencia numérica = 0.77855513624879024714

disp("--> Gauss-Seidel");
G=inv(D-C_L)*C_U;
[_, it, tasa]=GaussSeidel(A, b, 10^-7, sol);
disp(tasa, "Tasa de convergencia aproximada =");
// Tasa de convergencia numérica = 0.44444444206879080150

disp("--> Descenso Rápido");
G=inv(D-C_L)*C_U;
[_, it, tasa]=DescensoRapido(A, b, 10^-7, sol);
disp(tasa, "Tasa de convergencia aproximada =");
// Tasa de convergencia numérica = 0.71390754156507096884

disp("Se observa que las tasas numéricas son una buena aproximación de las tasas exactas. Si se imprimen las tasas iteración por iteración, se ve que esta esta oscilando alrededor de la tasa exacta.");



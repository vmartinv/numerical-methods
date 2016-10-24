funcprot(0);

/////////////////////////////////////////////////////////////////
disp("=============Ejercicio 1=============");
disp("Se implementa la funcion GaussSinPivoteo: Se convierte Ax=b en Ux=c, con U triangular superior.");
function [U,c]=GaussSinPivoteo(A,b)
  eps = 1e-10;
  n = length(b);
  for fila=1:n    
    if abs(A(fila,fila))<eps then printf("Error: se necesita hacer pivoteo"); abort; end
    // Eliminacion
    for i=fila+1:n
      z = A(i,fila) / A(fila,fila);
      A(i, fila) = 0;
      for j=fila+1:n
        A(i,j) = A(i,j) - z*A(fila,j);
      end
      b(i) = b(i) - z*b(fila);
    end  
  end
  c=b;
  U=A;
endfunction

disp("(+) Ejemplo del Kincaid, pagina 139:");
A=[6 -2 2 4; 12 -8  6 10; 3 -13 9 3; -6 4 1 -18];
b=[12; 34; 27; -38];
disp(A, "A = ");
disp(b, "b = ");
[U, c]=GaussSinPivoteo(A, b);
disp(U, "U = ");
disp(c, "c = ");


/////////////////////////////////////////////////////////////////
disp("=============Ejercicio 2=============");
disp("Se implementa la funcion GaussPivoteo: Se convierte Ax=b en Ux=c, con U triangular superior.");
function [U,c]=GaussPivoteo(A,b)
  eps = 1e-10;
  n = length(b);
  permut(n) = 0;
  for i=1:n permut(i) = i; end
  
  for fila=1:n
    // Pivoteo
    if abs(A(fila,fila))<eps then
      nfila = -1;
      for i=fila+1:n
        if abs(A(i, fila))>eps then nfila = i; break; end
      end
      if nfila == -1 then printf("Error: matriz singular!"); abort; end
      if nfila <> fila then
        for j=fila:n
          vaux = A(nfila,j);
          A(nfila,j) = A(fila,j);
          A(fila,j) = vaux;
        end
        vaux = b(nfila);
        b(nfila) = b(fila);
        b(fila) = vaux;
        vaux = permut(nfila);
        permut(nfila) = permut(fila);
        permut(fila) = vaux;
      end
    end
    
    // Eliminacion
    for i=fila+1:n
      z = A(i,fila) / A(fila,fila);
      A(i, fila) = 0;
      for j=fila+1:n
        A(i,j) = A(i,j) - z*A(fila,j);
      end
      b(i) = b(i) - z*b(fila);
    end
  end
  U=A;
  c=b;
endfunction

disp("(+) Ejemplo del Kincaid, pagina 143:");
A=[0 1; 1 1];
b=[1; 2];
disp(A, "A = ");
disp(b, "b = ");
[U, c] = GaussPivoteo(A, b);
disp(U, "U = ");
disp(c, "c = ");

disp("(+) Otro ejemplo:");
A=[0 1 0; 1 0 0; 0 0 1];
b=[1; 2; 3];
disp(A, "A = ");
disp(b, "b = ");
[U, c] = GaussPivoteo(A, b);
disp(U, "U = ");
disp(c, "c = ");

/////////////////////////////////////////////////////////////////
disp("=============Ejercicio 3=============");
disp("Se implementa la funcion Gauss: Se convierte Ax=b en Ux=c, con U triangular superior.");
function [U,c]=Gauss(A,b)
  eps = 1e-10;
  n = length(b);
  permut(n) = 0;
  for i=1:n permut(i) = i; end
  
  for fila=1:n
    // Pivoteo
    nfila = fila;
    nval = -1;
    for i=fila:n
      s=abs(A(i,fila));
      for j=fila+1:n
        s=max(s, abs(A(i, j)));
      end
      val = abs(A(i,fila)/s);
      if nval == -1 | val > nval then nfila = i;  nval = val; end
    end
    if abs(A(nfila,fila))<eps then printf("Error: matriz singular!"); abort; end
    if nfila <> fila then
      for j=fila:n
        vaux = A(nfila,j);
        A(nfila,j) = A(fila,j);
        A(fila,j) = vaux;
      end
      vaux = b(nfila);
      b(nfila) = b(fila);
      b(fila) = vaux;
      vaux = permut(nfila);
      permut(nfila) = permut(fila);
      permut(fila) = vaux;
    end
    
    // Eliminacion
    for i=fila+1:n
      z = A(i,fila) / A(fila,fila);
      A(i, fila) = 0;
      for j=fila+1:n
        A(i,j) = A(i,j) - z*A(fila,j);
      end
      b(i) = b(i) - z*b(fila);
    end
  end
  U=A;
  c=b;
endfunction


/////////////////////////////////////////////////////////////////
disp("=============Ejercicio 4=============");
disp("Se muestran a continuacion los errores cometidos respecto a la solución exacta.");
function x=Solve(A, b, gauss)
  [U, c] = gauss(A, b);
  n = length(c);
  x = zeros(n, 1);
  for i=n:-1:1
    x(i) = c(i);
    for j=i+1:n
      x(i) = x(i) - U(i, j)*x(j);
    end
    x(i) = x(i) / U(i, i);
  end
endfunction

function test(A, b, sol)
  disp(Solve(A, b, GaussPivoteo) - sol, "Gauss con pivoteo cualquiera:");
  disp(Solve(A, b, Gauss) - sol, "Gauss con criterio de eleccion de pivote:");
  disp(inv(A)*b - sol, "Usando matriz inversa:");
endfunction
disp("a)");
A=[10^-12 1;1 1];
b=[1; 2];
sol=[1/(1 - 10^-12); (1 - 2*10^-12)/(1 - 10^-12)];
test(A, b, sol);

disp("b)");
A=[4 5 -6;2 0 -7; -5 -8 0];
b=[-28; 29; -64];
sol=[5200/47; -2874/47; 1291/47];
test(A, b, sol);

disp("c)");
A=[1 2 -1 0 0 3 1; 1 2 2 1 -4 1 0; 0 1 -1 3 -3 0 0; 0 1 -1 2 1 1 0; 0 0 1 -2 1 0 1; 0 0 0 2 0 0 3; 0 0 0 1 1 -1 0];
b=[-2; -2; 2; 5; -7; -8; 2];
sol=[1; -1; 0; 2; 1; 1; -4];
test(A, b, sol);





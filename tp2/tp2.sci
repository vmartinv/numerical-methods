//LCC - Metodos Numericos
//Trabajo Practico 2
//Integrantes:
//    Matias Saper
//    Martin Villagra
//    Ciro Barbagelata

funcprot(0);

/////////////////////////////////////////////////////////////////
disp("=============Ejercicio 1=============");
disp("Se implementa la funcion GaussSinPivoteo: Se convierte Ax=b en Ux=c, con U triangular superior.");
function [U,c]=GaussSinPivoteo(A,b)
  eps = %eps;
  n = length(b);
  for fila=1:n    
    if abs(A(fila,fila))<eps then error("Se necesita hacer pivoteo"); end
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
disp("Se implementa la funcion GaussConPivoteo: Se convierte Ax=b en Ux=c, con U triangular superior.");
function [U,c]=GaussConPivoteo(A,b,eps)
  if ~exists("eps","local") then eps = %eps; end
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
      if nfila == -1 then error("Matriz singular!"); end
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
A=[10^-12 1; 1 1];
b=[1; 2];
disp(A, "A = ");
disp(b, "b = ");
[U, c] = GaussConPivoteo(A, b, 10^-10);
disp(U, "U = ");
disp(c, "c = ");

disp("(+) Otro ejemplo:");
A=[0 1 0; 1 0 0; 0 0 1];
b=[1; 2; 3];
disp(A, "A = ");
disp(b, "b = ");
[U, c] = GaussConPivoteo(A, b);
disp(U, "U = ");
disp(c, "c = ");

/////////////////////////////////////////////////////////////////
disp("=============Ejercicio 3=============");
disp("Se implementa la funcion GaussConPivoteoSelectivo: Se convierte Ax=b en Ux=c, con U triangular superior.");
function [U,c]=GaussConPivoteoSelectivo(A,b,eps)
  if ~exists("eps","local") then eps = %eps; end
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
    if abs(A(nfila,fila))<eps then error("Matriz singular!"); end
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
disp("Se muestran a continuacion los errores cometidos respecto a la solución exacta y su explicación.");
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
  try
    disp(abs(Solve(A, b, GaussSinPivoteo) - sol), "Usando GaussSinPivoteo:");
  catch
    disp("No es posible realizar el algoritmo sin pivoteo por quedarnos con el 0 como pivote.");
  end
  disp(abs(Solve(A, b, GaussConPivoteo) - sol), "Usando GaussConPivoteo:");
  disp(abs(Solve(A, b, GaussConPivoteoSelectivo) - sol), "Usando GaussConPivoteoSelectivo:");
  disp(abs(inv(A)*b - sol), "Usando matriz inversa:");
  disp(norm(A)*norm(inv(A)), "Numero de condicion de A:");
endfunction
 
disp("a)");
A=[10^-12 1;1 1];
b=[1; 2];
sol=[1/(1 - 10^-12); (1 - 2*10^-12)/(1 - 10^-12)];
test(A, b, sol);
// Gauss sin pivoteo:   
// 0.0000221  
// 0. 

// Gauss con pivoteo cualquiera:   
// 0.0000221
// 0.  

// Gauss con criterio de eleccion de pivote:   
// 0.  
// 0. 
 
// Usando matriz inversa:   
// 0.         
// 1.110D-16  

// Numero de condicion de A:   
// 2.618034  

disp("Tenemos un error considerable en el Gauss sin pivoteo, prodecemos a analizarlo. Resolviendo el sistema tenemos:");
disp("x2 = (2 - (10^-12)^-1) / (1 - (10^-12)^-1), x1 = (1-x2)*(10^-12)^-1)");
disp("El termino 2-(10^-12)^-1 sera calculado como -(10^-12)^-1, pues al hacer el shift de las mantisas los bits de 2 se eliminan. Analogamente el denominador 1-(10^-12)^-1 tambien sera computado como -(10^-12)^-1. Por lo tanto x2 = 1 y x1 sera computado como 0, lo cual introduce error en la solucion.");
disp("Este error no se presenta cuando el pivote se elige cuidadosamente, pues esta elección reduce los errores de las operaciones. También observamos que el error usando la matriz inversa no es grande debido al bajo numero de condición que presenta A.");


disp("b)");
A=[4 5 -6;2 0 -7; -5 -8 0];
b=[-28; 29; -64];
sol=[5200/47; -2874/47; 1291/47];
test(A, b, sol);
// Gauss sin pivoteo:   
// 10^(-13) *
//    0.1421085  
//    0.0710543  
//    0.0355271  

// Gauss con pivoteo cualquiera:   
// 10^(-13) *
//    0.1421085  
//    0.0710543  
//    0.0355271  

// Gauss con criterio de eleccion de pivote:   
// 10^(-13) *
//    0.1421085  
//    0.0710543  
//    0.0355271  

// Usando matriz inversa:   
// 10^(-13) *
//    0.2842171  
//    0.1421085  
//    0.1065814  

// Numero de condicion de A:   
// 26.094814

disp("Todos los algoritmos producen el mismo error, es decir que en este caso realizar un pivoteo solo empeora el tiempo de ejecución sin obtener mejores resultados.");
disp("Tenemos un problema al usar la matrix inversa, donde los errores casi se duplican. Veamos que el numero de condicion para la matriz A es alto, por lo que el error obtenido es esperable.");


disp("c)");
A=[1 2 -1 0 0 3 1; 1 2 2 1 -4 1 0; 0 1 -1 3 -3 0 0; 0 1 -1 2 1 1 0; 0 0 1 -2 1 0 1; 0 0 0 2 0 0 3; 0 0 0 1 1 -1 0];
b=[-2; -2; 2; 5; -7; -8; 2];
sol=[1; -1; 0; 2; 1; 1; -4];
test(A, b, sol);
//  No es posible realizar el algoritmo sin pivoteo por quedarnos con el 0 como pivote.

// Gauss con pivoteo cualquiera:   
// 10^(-14) *
//    0.3552714  
//    0.1776357  
//    0.         
//    0.0888178  
//    0.0222045  
//    0.         
//    0.         

// Gauss con criterio de eleccion de pivote:   
// 0.         
// 1.221D-15  
// 2.961D-16  
// 4.441D-16  
// 2.220D-16  
// 5.551D-16  
// 8.882D-16  

// Usando matriz inversa:   
// 10^(-14) *
//    0.1776357  
//    0.0888178  
//    0.0777156  
//    0.1110223  
//    0.         
//    0.1110223  
//    0.0888178

// Numero de condicion de A:   
// 18.685037

disp("Aquí también se aprecia error en el calculo de la inversa debido a su valor de condición. A su vez tenemos un mejor resultado cuando utilizamos Gauss. Podría parecer que GaussConPivoteoSelectivo es peor pero al analizar la norma del error, notamos que utilizando GaussConPivoteoSelectivo la norma del error se reduce a más de la mitad con respecto a GaussConPivoteo. Esto se debe a que el criterio reduce la propagación del error en el algoritmo.");

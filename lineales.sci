// Algoritmo de Gauss modificado con pivoteo parcial.
// Obtiene también la descomposición PA=LU ;) ;) (con L(i,i)=1)
function [x, P, L, U]=Gauss(A,b)
  n = size(A, 'r');
  if ~exists("b","local") then
    for i=1:n b(i)=0; end
  end
  L=zeros(n,n);
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
    if abs(A(nfila,fila))<%eps then error("Matriz singular!"); end
    if nfila <> fila then
      for j=fila:n
        vaux = A(nfila,j);
        A(nfila,j) = A(fila,j);
        A(fila,j) = vaux;
      end
      for j=1:fila
        vaux = L(nfila,j);
        L(nfila,j) = L(fila,j);
        L(fila,j) = vaux;
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
      L(i,fila) = z;
      for j=fila+1:n
        A(i,j) = A(i,j) - z*A(fila,j);
      end
      b(i) = b(i) - z*b(fila);
    end
  end
  P = zeros(n,n);
  for i=1:n 
    P(permut(i),i)=1;
  end
  for i=1:n L(i,i)=1; end
  U=A;
  x=b;
endfunction

// Resuelve el sistema L*x = b,
// siendo L triangular inferior con diagonal no nula (inversible)
function x=SustitucionProgresiva(L,b)
  n = length(b);
  for i=1:n
    x(i) = b(i);
    for j=1:i-1
      x(i) = x(i) - L(i,j)*x(j);
    end
    x(i) = x(i) / L(i,i);
  end
endfunction

// Resuelve el sistema U*x = b,
// siendo U triangular superior con diagonal no nula (inversible)
function x=SustitucionRegresiva(U,b)
  n = length(b);
  for i=n:-1:1
    x(i) = b(i);
    for j=n:-1:i+1
      x(i) = x(i) - U(i,j)*x(j);
    end
    x(i) = x(i) / U(i,i);
  end
endfunction

// Resuelve Ax=b, dado PA=LU
function x=SolvePALU(b, P, L, U)
  z = SustitucionProgresiva(L, P*x);
  x = SustitucionRegresiva(U, z);
endfunction

// Resuelve el sistema A*x = b, siendo A tridiagonal
function x=Tridiagonal(A, b)
  n = length(b);
  for i=1:n d(i)=A(i, i); end
  for i=1:n-1 a(i)=A(i+1, i); end
  for i=1:n-1 c(i)=A(i, i+1); end
  for i=2:n
    d(i) = d(i) - (a(i-1) / d(i-1))*c(i-1);
    b(i) = b(i) - (a(i-1) / d(i-1))*b(i-1);
  end
  x(n) = b(n)/d(n);
  for i=n-1:-1:1
    x(i) = (b(i) - c(i)*x(i+1)) / d(i);
  end
endfunction

// Obtiene la factorización de Doolittle (A=LU, L(i, i)=1)
function [L, U]=Doolittle(A)
  n = size(A, 'r');
  for k=1:n
    L(k, k) = 1;
    for j=k:n
      sum = 0;
      for s=1:k-1
        sum = sum + L(k, s)*U(s, j);
      end
      U(k, j) = A(k, j) - sum
    end
    for i=k+1:n
      sum = 0;
      for s=1:k-1
        sum = sum + L(i, s)*U(s, k);
      end
      L(i, k) = (A(i, k) - sum) / U(k, k);
    end
  end
endfunction

// Obtiene la factorización de Crout (A=LU, U(i, i)=1)
function [L, U]=Crout(A)
  n = size(A, 'r');
  for k=1:n
    U(k, k) = 1;
    for i=k:n
      sum = 0;
      for s=1:k-1
        sum = sum + L(i, s)*U(s, k);
      end
      L(i, k) = A(i, k) - sum
    end
    for j=k+1:n
      sum = 0;
      for s=1:k-1
        sum = sum + L(k, s)*U(s, j);
      end
      U(k, j) = (A(k, j) - sum) / L(k, k);
    end
  end
endfunction

// Descomposicion de Cholesky, si A es real, simétrica y def. +. (A=L L^T)
// No necesita pivoteo por que sus elementos estan acotados por los de A.
function L=Cholesky(A)
  n = size(A, 'r');
  for k=1:n
    sum = 0;
    for s=1:k-1
      sum = sum + L(k, s)^2;
    end
    L(k, k) = sqrt(A(k, k) - sum);
    for i=k+1:n
      sum = 0;
      for s=1:k-1
        sum = sum + L(i, s)*L(k, s);
      end
      L(i, k) = (A(i, k) - sum) / L(k, k);
    end
  end
endfunction

// Descompone A = D - C_L - C_U con D diagonal,
// C_L negativo de la parte triangular inferior
// y C_U negativo de la parte triangular superior,
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

// Resuelve Ax=b usando Q=I, G=I-A
// Converge si es diagonalmente dominante con a_ii=1 o ||I-Q^-1 A||<1
function x=Richardson(A, b, it, x) //No testeado!
  n = length(b);
  if ~exists("x","local") then
    for i=1:n x(i)=0; end
  end
  for k=1:it
    for i=1:n
      sum = 0;
      for j=1:n
        sum = sum + A(i, j)*x(j);
      end
      r(i) = b(i) - sum;
    end
    for i=1:n
      x(i) = x(i) + r(i);
    end
  end
endfunction

// Resuelve Ax=b usando Q=D, G=D^-1(C_L+C_U)
// Converge si es diagonalmente dominante o ||I-Q^-1 A||<1
function x=Jacobi(A, b, it, x)
  n = length(b);
  if ~exists("x","local") then
    for i=1:n x(i)=0; end
  end
  for k=1:it
    for i=1:n
      sum = 0;
      for j=1:n
        if j<>i then
          sum = sum + A(i, j)*x(j);
        end
      end
      u(i) = (b(i) - sum) / A(i, i);
    end
    for i=1:n
      x(i) = u(i);
    end
  end
endfunction

// Resuelve Ax=b usando Q=D-C_L, G=(D-C_L)^-1 C_U
// Converge si es diagonalmente dominante o ||I-Q^-1 A|| < 1
function x=GaussSeidel(A, b, it, x)
  n = length(b);
  if ~exists("x","local") then
    for i=1:n x(i)=0; end
  end
  for k=1:it
    for i=1:n
      sum = 0;
      for j=1:n
        if j<>i then
          sum = sum + A(i, j)*x(j);
        end
      end
      x(i) = (b(i) - sum) / A(i, i);
    end
  end
endfunction

function C=dot(A, B)
  C=A'*B;
endfunction

// Resuelve Ax=b, minimizando q(x)=<x,Ax>-2<x,b>,
// moviendose en direccion opuesta a su gradiente (2(Ax-b))
function x=DescensoRapido(A, b, it, x)
  n = length(b);
  if ~exists("x","local") then
    for i=1:n x(i)=0; end
  end
  for k=1:it
    v = b - A*x;
    div = dot(v, A*v);
    if abs(div)<%eps then break; end
    t = dot(v, v) / div;
    x = x + t*v;
 end
endfunction

// Retorna el autovalor de mayor valor absoluto (r),
// junto con un autovector asociado (x).
function [x, r]=MetodoPotencia(A, it, x)
  n = size(A, 'r');
  if ~exists("x","local") then
    for i=1:n x(i)=1; end
  end
  for k=1:it
    y = A*x;
    r = max(y)/max(x); // Tambien se puede probar con r = y(1)/x(1)
    x = y/norm(y);
 end
endfunction

// Retorna el autovalor de mayor valor absoluto (r),
// junto con un autovector asociado (x).
// Funciona aplicando el método de la potencia sobre A^-1.
function [x, r]=MetodoPotenciaInversa(A, it, x)
  n = size(A, 'r');
  if ~exists("x","local") then
    for i=1:n x(i)=1; end
  end
  [_, P, L, U]=Gauss(A);
  for k=1:it
    y = SolvePALU(x, P, L, U);
    r = max(x)/max(y); // Tambien se puede probar con r = x(1)/y(1)
    x = y/norm(y);
 end
endfunction



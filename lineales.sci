// Algoritmo de Gauss modificado con pivoteo parcial.
// Obtiene también la descomposición PLU ;) ;) (con L(i,i)=1)
function [U, x, L, P]=Gauss(A,b)
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

//Obtiene la factorización de Doolittle (A=LU), l_ii=1
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

//Obtiene la descomposicion de Cholesky, si A es real, simétrica y definida positiva. (A=L L^T)
//Esta factorizacion no necesita pivoteo por que sus elementos estan acotados por los de A.
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

// Resuelve el sistema U*x = b, siendo U triangular superior con diagonal no nula (inversible)
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

// Resuelve el sistema L*x = b, siendo L triangular inferior con diagonal no nula (inversible)
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

//Resuelve Ax=b usando Q=I, G=I-A
//x^(k) = (I - A) x^(k-1) + b
function x=Richardson(A, b, it, x) //No testeada
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

//Resuelve Ax=b usando Q=diagonal de A, G=D^-1(C_L+C_U)
//Dx^(k) = (C_L + C_U) x^(k-1) + b
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

//Resuelve Ax=b usando Q=parte triangular inferior de A, G=(D-C_L)^-1 C_U
//(D-C_L)x^(k) = C_U x^(k-1) + b
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


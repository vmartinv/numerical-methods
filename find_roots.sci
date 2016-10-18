//// Metodo de la biseccion
function c = biseccion(f,a,b,its,eps,delta)
    if ~exists("eps","local") then eps = %eps; end
    if ~exists("delta","local") then delta = %eps; end
    if sign(f(a)) == sign(f(b)) then printf("Bad preconditions"); abort end
    for i=1:its do
        c = a+(b-a)/2;
        if abs(b-a) < delta | abs(f(c)) < eps then break; end
        if sign(f(a)) == sign(f(c)) then a=c
        else b=c;
        end
    end
endfunction

//// Devuelve el numero de iteraciones necesarias para alcanzar tal error
function its = biseccionError(a,b,err)
    its = ceil(max(0,log2((b-a)/abs(a)/err)-1));
endfunction

//// Metodo de Newton
// Uso con derivada numerica para f(x)=x^3+e^x
// deff('r=f(x)','r=x^3+%e^x')
// deff('r=df(x)','r=numderivative(f,x,0.0001)')
// newton(f,df,x0,100)
// Sirve también para multivariable donde df(c) es la matriz jacobiana de f en c
// Uso con derivada numerica para 0=x^2+x*y^3-9 0=3*x^2*y-4-y^3
// deff('r=f(x)','r=[x(1)^2+x(1)*x(2)^3-9;3*x(1)^2*x(2)-4-x(2)^3]')
// deff('r=jf(x)','r=numderivative(f,x,0.0001)')
function xi=newton(f,df,x0,its,eps)
    if ~exists("eps","local") then eps = %eps; end
    xi = x0;
    for i=1:its
        xi = xi - inv(df(xi))*f(xi);
        if abs(x0-xi)<eps then break; end
        x0 = xi
    end
endfunction

//// Metodo de la secante
function xn=secante(f,x0,x1,its,eps,delta)
    if ~exists("eps","local") then eps = %eps; end
    if ~exists("delta","local") then delta = %eps; end
    xn_1 = x0; fxn_1 = f(xn_1);
    xn = x1; fxn = f(xn);
    for i=1:its
        if abs(fxn_1) < abs(fxn) then //swap xn y xn_1
            tmp = xn
            xn = xn_1
            xn_1 = tmp
            tmp = fxn 
            fxn = fxn_1
            fxn_1 = tmp
        end
        if abs(fxn-fxn_1)<eps then printf("Warning: funcion monotona"); break; end
        xnext = xn - fxn * (xn-xn_1)/(fxn-fxn_1);
        xn_1 = xn; fxn_1 = f(xn_1);
        xn = xnext; fxn = f(xn);
        if abs(fxn) < eps | abs(xn-xn_1) < delta then break; end
    end
endfunction

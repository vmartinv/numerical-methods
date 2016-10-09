//// Metodo de la biseccion
function r = biseccion(f,a,b,its)
    if sign(f(a)) == sign(f(b)) then printf("Bad preconditions"); abort end
    for i=1:its do
        c = a+(b-a)/2;
        if sign(f(a)) == sign(f(c)) then a=c
        else b=c;
        end
    end
    r = a+(b-a)/2;
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
function xi=newton(f,df,x0,its,err)
    if ~exists("err","local") then
        err = 0
    end
    xi = x0;
    for i=1:its
        xi = xi - inv(df(xi))*f(xi);
        if abs(x0-xi)<err then break;
        end
        x0 = xi
    end
endfunction

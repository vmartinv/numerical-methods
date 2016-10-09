exec('./find_roots.sci',-1)

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 1.2==")
deff('r=f(x)','r=x^3+1')
disp("x="+string(biseccion(f, -1, 10, 40)))

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 1.8==")
disp("Se necesitan "+string(biseccionError(-1,10,0.00000001))+" iteraciones.")

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 1.9==")
deff('r=f(x)','r=x^-1 - tan(x)')
try
    biseccion(f, 0, %pi/2, 20)
catch
    disp("a) f no es continua en 0.")
end

deff('r=f(x)','r=x^-1 - 2^x')
try
    biseccion(f, 0, 1, 20)
catch
    disp("b) f no es continua en 0.")
end

deff('r=f(x)','r=2^-x + %e^x + 2*cos(x) -6')
disp("c) x="+string(biseccion(f, 1, 3, 200)))

deff('r=f(x)','r=2*x^3 - 9*x^2 + 18*x -2')
disp("d) f no es continua en "+string(biseccion(f, 0, 4, 200)))

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 1.12==")
disp("a) Verdadero")
disp("b) Verdadero")
disp("c) Falso, ejemplo: a_n=2, b_n=4, r=3+EPS")
disp("d) Verdadero")
disp("e) Falso, a_n=2, b_n=4, r=4-EPS")
disp("f) False, mismo ejemplo que el c).")

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 1.16==")
a=0.1
b=0.1
disp("(a+b)/2="+string((a+b)/2), "a+0.5*(b-a)="+string(a+0.5*(b-a)), "b="+string(b), "a="+string(a))
disp("TODO")

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 1.19==")
deff('r=f(x)','r=x - tan(x)')
disp("x="+string(biseccion(f, 1, 2, 200)))

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 2.1==")
deff('r=f(x)','r=tan(x)')
deff('r=df(x)','r=numderivative(f,x,0.0001)')
disp("x="+string(newton(f, df, 4.5, 10)))
disp("x="+string(newton(f, df, 7.7, 10)))

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 2.10==")
deff('r=f(x)','r=x^3+3*x-5*x^2-7')
deff('r=df(x)','r=numderivative(f,x,0.0001)')
disp("x="+string(newton(f, df, 5, 10)))

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 2.18==")
disp("TODO")

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 2.30==")
deff('r=f(x)','r=x^3+94*x^2-389*x+294')
deff('r=df(x)','r=numderivative(f,x,0.0001)')
disp("x="+string(newton(f, df, 2, 10)))
x=[-100:0.01:100]
//plot(x, f(x))
disp("La derivada en 2 tiene un valor muy bajo (descomentar el plot). Por lo que la correccion va a ser exageradamente grande.")

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 2.33==")
deff('r=f(x)','r=[1+(x(1))^2+(x(2))^2+(%e^x(1))*cos(x(2)); 2*x(1)*x(2)+%e^x(1)*sin(x(2))]')
deff('r=jf(x)','r=numderivative(f,x,0.0001)')
x=[-10:0.01:10]
xx = linspace(-10,10,40);
yy = linspace(-10,10,40);

function z=f2(x, y)
    a=[x;y];
    disp("x="+string(x))
    z=f(a)(1)
    disp(z)
endfunction
//deff('z=f2(x, y)', 'z=(f([x;y]))(3)')
zz = feval(xx, yy, f2)';
surf(xx,yy,zz)

disp("Notamos que el metodo tarda muchas iteraciones en converger")

disp("Con 10 iteraciones [x, y]=")
x=newton(f, jf, [-1; 4], 10)
disp(x)
disp(f(x))
disp("Con 35 iteraciones [x, y]=")
x=newton(f, jf, [-1; 4], 35)
disp(x)
disp(f(x))
disp("Con 40 iteraciones [x, y]=")
x=newton(f, jf, [-1; 4], 40)
disp(x)
disp(f(x))

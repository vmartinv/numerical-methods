funcprot(0)
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
deff('r=f(x)','r=[1+(x(1))^2+(x(2))^2+(%e^x(1))*cos(x(2)); 2*x(1, 1)*x(2)+%e^x(1)*sin(x(2))]')
deff('r=jf(x)','r=numderivative(f,x,0.0001)')
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

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 3.a==")
function r=f(x)
    k1=x(1); k2=x(2); k3=x(3);
    deff('r=ec(p, r)', 'r = k1*%e^(k2*r) + k3*r - p')
    r=[ec(10, 1);ec(12, 2);ec(15, 3)];
endfunction
deff('r=jf(x)','r=numderivative(f,x,0.0001)')
rand('seed', 2342345)
while %T
    try
        x=newton(f, jf, [rand(); rand(); rand()], 50);
        if find(isnan(x)) then continue; else break; end
    catch
        continue
    end
end
k1=x(1); k2=x(2); k3=x(3);
disp("k1="+string(k1)+", k2="+string(k2)+", k3="+string(k3))

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 3.b==")
deff('d=f(r)', 'd = k1*%e^(k2*r) + k3*r - 500')
r=biseccion(f, 1, 500, biseccionError(1, 500, 10^-12))
disp('r>='+string(r)+' pulgadas')

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 4.3==")
deff('r=f(x)', 'r=sin(x/2)-1')
x=secante(f, 3, 3.18, 15)
disp('a) x='+string(x)+", f(x)="+string(f(x)))

deff('r=f(x)', 'r=%e^x - tan(x)')
x=secante(f, -4, 2, 10)
disp('b) x='+string(x)+", f(x)="+string(f(x)))

deff('r=f(x)', 'r=x^3-12*x^2+3*x+1')
//x=[0:0.01:50]; plot(x, f(x))
x=secante(f, -1, -2, 10)
disp('c) x='+string(x)+", f(x)="+string(f(x)))

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 4.7==")
deff('r=f(x)', 'r=2-(x-1)/2')
x=secante(f, 1, 2, 1)
disp('x_2='+string(x))

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 4.8.a)==")
deff('r=f(x)', 'r=atan(x)-2*x/(1+x^2)')
x=biseccion(f, -2, -1, 200, %eps/2)
disp('Met. biseccion: x='+string(x)+', f(x)='+string(f(x)))

deff('r=df(x)','r=numderivative(f,x,0.0001)')
x=newton(f, df, 0.5, 10)
disp('Met. newton: x='+string(x)+', f(x)='+string(f(x)))

x=secante(f, 0.5, 1, 10)
disp('Met. secante: x='+string(x)+', f(x)='+string(f(x)))

x=biseccion(f, -0.5, 0.5, 200)
disp('Met. biseccion: x='+string(x)+', f(x)='+string(f(x)))

deff('r=df(x)','r=numderivative(f,x,0.0001)')
x=newton(f, df, -0.1, 50, 0)
disp('Met. newton: x='+string(x)+', f(x)='+string(f(x)))

x=secante(f, -0.5, -0.1, 10, %eps)
disp('Met. secante: x='+string(x)+', f(x)='+string(f(x)))

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 4.8.b)==")
deff('r=f(x)', 'r=atan(x)')
deff('r=df(x)','r=numderivative(f,x,0.0001)')
x=newton(f, df, 0, 50, 0)
disp('x='+string(x)+', f(x)='+string(f(x)))

/////////////////////////////////////////////////////////////////
disp("==Ejercicio 4.11==")
disp("La formula es equivalente trivialmente. Es peor porque usa mas multiplicaciones. ")

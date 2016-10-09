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

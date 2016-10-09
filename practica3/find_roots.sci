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

function r = biseccion(f,a,b,its)
    if sign(f(a)) == sign(f(b)) then printf("Bad preconditions"); return; end
    for i=1:its do
        c = a+(b-a)/2;
        if sign(f(a)) == sign(f(c)) then a=c
        else b=c;
        end
    end
    r = a+(b-a)/2;
endfunction

function r=biseccionError(a,b,err)
    r=ceil(max(0,log2((b-a)/abs(a)/err)-1));
endfunction


funcprot(0);

/////////////////////////////////////////////////////////////////
disp("Ejercicio 2");
e12=0.000006144212353328210;

function err=etaylor(n, x)
    fact=1;
    r=0;
    err=[];
    for i=0:n
        r=r+x^i/fact;
        fact=fact * (i+1);
        if i>=1 then 
            err(i)(1) = i;
            err(i)(2) = abs(e12-r);
        end
    end
endfunction

disp("Utilizando el desarrollo de Taylor, se obtienen los siguientes errores:");
disp("   n       |Error");
disp(etaylor(100, -12));
//~ n       |Error
//~ 1.      11.000006  
//~ 2.      60.999994  
//~ 3.      227.00001  
//~ 4.      636.99999  
//~ 5.      1436.6     
//~ 6.      2710.6     
//~ 7.      4398.8857  
//~ 8.      6265.3429  
//~ 9.      7953.6286  
//~ 10.     9109.1371  
//~ 11.     9504.7891  
//~ 12.     9109.1371  
//~ 13.     8072.9486  
//~ 14.     6654.5535  
//~ 15.     5127.4482  
//~ 16.     3709.053   
//~ 17.     2528.4772  
//~ 18.     1629.8763  
//~ 19.     996.45226  
//~ 20.     579.34486  
//~ 21.     321.11064  
//~ 22.     170.04691  
//~ 23.     86.209204  
//~ 24.     41.918851  
//~ 25.     19.582615  
//~ 26.     8.8026771  
//~ 27.     3.8130084  
//~ 28.     1.5937139  
//~ 29.     0.6435505  
//~ 30.     0.2513553  
//~ 31.     0.0950599  
//~ 32.     0.0348458  
//~ 33.     0.0123926  
//~ 34.     0.0042798  
//~ 35.     0.0014365  
//~ 36.     0.0004690  
//~ 37.     0.0001490  
//~ 38.     0.0000461  
//~ 39.     0.0000139  
//~ 40.     0.0000041  
//~ 41.     0.0000012  
//~ 42.     0.0000003  
//~ 43.     9.042D-08  
//~ 44.     2.423D-08  
//~ 45.     6.348D-09  
//~ 46.     1.628D-09  
//~ 47.     4.082D-10  
//~ 48.     1.008D-10  
//~ 49.     2.383D-11  
//~ 50.     6.093D-12  
//~ 51.     9.478D-13  
//~ 52.     6.769D-13  
//~ 53.     3.091D-13  
//~ 54.     3.908D-13  
//~ 55.     3.730D-13  
//~ 56.     3.768D-13  
//~ 57.     3.760D-13  
//~ 58.     3.762D-13  
//~ 59.     3.761D-13  
//~ 60.     3.761D-13  
//~ 61.     3.761D-13  
//~ 62.     3.761D-13  
//~ 63.     3.761D-13  
//~ 64.     3.761D-13  
//~ 65.     3.761D-13  
//~ 66.     3.761D-13  
//~ 67.     3.761D-13  
//~ 68.     3.761D-13  
//~ 69.     3.761D-13  
//~ 70.     3.761D-13  
//~ 71.     3.761D-13  
//~ 72.     3.761D-13  
//~ 73.     3.761D-13  
//~ 74.     3.761D-13  
//~ 75.     3.761D-13  
//~ 76.     3.761D-13  
//~ 77.     3.761D-13  
//~ 78.     3.761D-13  
//~ 79.     3.761D-13  
//~ 80.     3.761D-13  
//~ 81.     3.761D-13  
//~ 82.     3.761D-13  
//~ 83.     3.761D-13  
//~ 84.     3.761D-13  
//~ 85.     3.761D-13  
//~ 86.     3.761D-13  
//~ 87.     3.761D-13  
//~ 88.     3.761D-13  
//~ 89.     3.761D-13  
//~ 90.     3.761D-13  
//~ 91.     3.761D-13  
//~ 92.     3.761D-13  
//~ 93.     3.761D-13  
//~ 94.     3.761D-13  
//~ 95.     3.761D-13  
//~ 96.     3.761D-13  
//~ 97.     3.761D-13  
//~ 98.     3.761D-13  
//~ 99.     3.761D-13  
//~ 100.    3.761D-13  

disp("Vemos que en las primeras iteraciones hay una gran variacion en el resultado.");
disp("El error se puede observar mirando el error del polinomio de Taylor de grado n:");
disp("e^eps * (-12)^n / (n+1)! con eps entre 0 y -12");
disp("Para valores pequeños de n, el termino (-12)^n va a ser muy grande, mientras que e^eps / (n+1)! va a ser relativamente pequeño.");


// Falta dar un metodo alternativo (plantear la ec. y aplicar logaritmo a ambos lados.) 


/////////////////////////////////////////////////////////////////
disp("Ejercicio 3");

disp("Ejercicio 3.a");
deff('r=f(x)','r=x^3-log(1+2*x)');
disp("Se muestra el grafico de la funcion, donde se ve una raiz cercana al 1.");
plot([1:0.001:2],f);


disp("Ejercicio 3.b");
disp("Se puede usar el metodo de Steffensen, el cual viene dado por:")
disp("x_{n+1} = x_n + f(x_n)/g(x_n)")
disp(", donde g(x) = (f(x_n+f(x_n)) - f(x_n))/f(x_n)");
disp("La funcion g es la pendiente de la funcion entre (x_n, f(x_n)) y el punto auxiliar (x_n+h, f(x_n+h), con h=f(x_n).");
disp("Cuando se provee un valor lo suficientemente cercano a la raiz, h va a ser pequeño y g(x_n) es una buena aproximacion de f''(x_n).");
disp("En otras palabras, cuando el valor es elegido correctamente el método se comporta como el de Newton, el cual converge por encontrarse el x inicial cerca de la raíz deseada.");

function xi=steff(f,x0,its,eps)
    if ~exists("eps","local") then eps = %eps; end
    if abs(f(x0))<eps then
        xi=x0;
        return;
    end
    xi_1 = x0;
    for i=1:its
        xi = xi_1 - f(xi_1)^2/(f(xi_1+f(xi_1)) - f(xi_1));
        if abs(xi_1-xi)<eps | abs(f(xi))<eps then break; end
        xi_1 = xi;
    end
endfunction

r=steff(f,1,200);
// r  =    1.0400253

disp("Ejercicio 3.c");
function tab=tabla(f,x0,r)
    tab=[];
    for j=4:20
        tab(j-3)(1)=10^-j
        for i=1:100
            if abs(steff(f,x0,i) - r) <= 10^-j then
                tab(j-3)(2)=i;
                break;
             end
        end
   end
endfunction
disp("Error           |Iteraciones");
disp(tabla(f,1,r));

//~ Error           |Iteraciones   
//~ 0.0001       3.  
//~ 0.00001      3.  
//~ 0.000001     3.  
//~ 1.000D-07    4.  
//~ 1.000D-08    4.  
//~ 1.000D-09    4.  
//~ 1.000D-10    4.  
//~ 1.000D-11    4.  
//~ 1.000D-12    5.  
//~ 1.000D-13    5.  
//~ 1.000D-14    5.  
//~ 1.000D-15    5.  
//~ 1.000D-16    5.  
//~ 1.000D-17    5.  
//~ 1.000D-18    5.  
//~ 1.000D-19    5.  
//~ 1.000D-20    5. 

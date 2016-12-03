/////////////////////////////////////////////////////////////////
disp("=============Ejercicio 2.(a)=============");
// La implementación de Newton utilizando las diferencias divididas
// adaptado del Kincaid, página 309 (Cap. 6.2)
// Se agrego +1 en todas las indexaciones.
function d=DifDivididasNewton(xi,yi)
    n = length(xi)-1;
    for i=0:n
        d(i+1) = yi(i+1);
    end
    for j=1:n
      for i=n:-1:j
        d(i+1) = (d(i+1) - d(i)) / (xi(i+1) - xi(i-j+1));
      end
    end
endfunction

// Sirve para evaluar los coeficientes que devuelve Newton
function s=EvalNewton(xi,d,x)
    n=length(xi)
    p = 1; s = 0;
    for i=1:n
        s = s + d(i)*p;
        p = p * (x-xi(i));
    end
endfunction


xi=[-1 1 2 4 10 18 25];
yi=[10 11.5 13.3 24.7 101 64.2 159];
d=DifDivididasNewton(xi, yi);
for i=1:length(xi)
  E(i) = abs((yi(i)-EvalNewton(xi,d,xi(i)))/yi(i));
end
disp(max(E), "Máximo error relativo de la interpolación sobre los puntos: ");
// Máximo error relativo de la interpolación sobre los puntos: 5.720D-15

// Devuelve una string representando el polinomio en su forma de Newton
function s=NewtonString(xi,d)
    n=length(xi)
    p = ""; s = "";
    for i=1:n
        if length(s)>0 then s = s + " + "; end
        s = s + string(d(i));
        if length(p)>0 then s = s + "*" + p; end
        if length(p)>0 then p = p + "*"; end
        p = p + "(x"
        if xi(i)<0 then p = p + "+";
        elseif xi(i)>0 then p = p + "-";
        end
        if xi(i)<>0 then p = p + string(abs(xi(i))); end
        p = p + ")";
    end
endfunction
disp(NewtonString(xi,d), "f(x) ~= ");
// f(x) ~= 10 + 0.75*(x+1) + 0.35*(x+1)*(x-1) + 0.19*(x+1)*(x-1)*(x-2) +
// -0.0215446*(x+1)*(x-1)*(x-2)*(x-4) + 0.0008704*(x+1)*(x-1)*(x-2)*(x-4)*(x-10)
// + -0.0000081*(x+1)*(x-1)*(x-2)*(x-4)*(x-10)*(x-18)

/////////////////////////////////////////////////////////////////
disp("=============Ejercicio 2.(b)=============");
disp(EvalNewton(xi,d,0), "f(0) ~= ");
// f(0) ~= 11.033688

/////////////////////////////////////////////////////////////////
disp("=============Ejercicio 3.(c)=============");
disp("Se utilizan los censos de la ciudad de Rosario.");
disp("No se disponen de 10 pares debido a que sólo se encontraron 7 censos. Fuera de los censos existen estimaciones año a año que utilizan otra función de aproximación distinta a la aquí usada.");
xi=[1947 1960 1970 1980 1991 2001 2010];
yi=[484021 626845 750455 794127 894645 908163 948312];
// Fuentes utilizadas:
// América Latina: Urbanización y Evolución de la Población Urbana 1950-2000, CEPAL/CELADE.
// Población de la ciudad de Rosario, Dirección General de Estadística. Elaboración propia en base a datos INDEC-IPEC y Dirección Gral. de Topografía y Catastro.

function [c, a, f]=MinCuadsExpFit(xi, yi)
  n = length(xi);
  yi = log(yi);
  coci = n*sum(xi.^2) - sum(xi)^2
  a = (n*sum(xi.*yi) - sum(xi)*sum(yi)) / coci;
  b = (sum(xi.^2)*sum(yi) - sum(xi)*sum(xi.*yi)) / coci;
  c = %e^b;
  deff('[y]=f(x)','y=c*%e^(a*x)');
endfunction

[c, a, f] = MinCuadsExpFit(xi, yi);
E = [];
for i=1:length(xi)
  E(i) = abs((yi(i)-f(xi(i)))/yi(i));
end
disp(max(E), "Máximo error relativo de la aproximación sobre los puntos: ");
// Máximo error relativo de la aproximación sobre los puntos: 0.1160555.
disp(round(f(2020)), "Población estimada para el 2020: ");
// Población estimada para el 2020: 1135190.

disp("Este resultado es exagerado. Según la IPEC, para el 2020 se espera contar con apenas un millón de habitantes (por lo que el censo debería tener un recuento menor a eso). Este exageración posiblemente se deba a que Rosario no está creciendo como ciudad, puesto que cada vez más y más personas van a vivir a ciudades del Gran Rosario, como Funes (que experimento un crecimiento de casi el doble de su población en la última década). Este comportamiento no ocurría durante los primeros años de la ciudad, por lo que es difícil reflejar este cambio con una aproximación exponencial. De hecho si quitamos los primeros años al hacer el método la predicción se vuelve más razonable.");
// Fuente: Población estimada al 1° de julio de cada año calendario, según departamento y localidad. Provincia de Santa Fe. Años 2010-2025, IPEC.

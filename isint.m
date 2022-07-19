function f = isint(x)
epsi = 1e-6;
f = abs(x-round(x)) < epsi;

function y = morlet_m(f,t,width)
sf = f/width;
st = 1/(2*pi*sf);
%A = 1/sqrt(2*pi*st^2);
A=1/sqrt(st*sqrt(pi));
y = A*exp(-t.^2/(2*st^2)).*exp(i*2*pi*f.*t);
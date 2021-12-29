%µÍÍ¨ÂË²¨Æ÷
function [t,st]=lpf(f,sf,B)
df = f(2)-f(1);
T = 1/df;
hf = zeros(1,length(f));%È«Áã¾ØÕó
bf = [-floor( B/df ): floor( B/df )] + floor( length(f)/2 );
hf(bf)=1;
yf=hf.*sf;
[t,st]=F2T(f,yf);
st = real(st);
function sse = sseval_catchbond(x,fdata,tdata)
A = x(1);
b = x(2);
C = x(3);
d = x(4);

sse = sum((tdata - 1./(A*exp(-b*fdata)+C*exp(d*fdata))).^2);
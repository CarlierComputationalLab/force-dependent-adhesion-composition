% function sse = sseval_slip(x,tData,concData)
function sse = sseval_slip(x,fdata,rdata)
k_off_UL = x(1);
k = x(2);
% C = x(3);
% d = x(4);

sse = sum((rdata - (k_off_UL*exp(k*fdata))).^2);
% sse = sum((rdata - (log(k_off_UL)+(k*fdata))).^2);

% k_off = x(1);
% sse = sum((concData-(exp(-k_off*tData))).^2);
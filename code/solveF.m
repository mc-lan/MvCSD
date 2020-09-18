function F = solveF(W,Yl)
nl = size(Yl,1);
class = length(unique(Yl));
YYl = zeros(nl,class);
for i=1:nl
    c = Yl(i);
    YYl(i,c) = 1;
end
W = (abs(W)+abs(W'))/2;
L=diag(sum(W,2))-W;
Luu = L(nl+1:end,nl+1:end);
Lul = L(nl+1:end,1:nl);
Fu = -(Luu+eye(size(Luu,1))*0.001)\(Lul*YYl);
F = [YYl; Fu];
end
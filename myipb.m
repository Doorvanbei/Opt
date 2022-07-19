function [x,v,f] = myipb(c,iind,Al,bl,A,b,lb,ub)
% MILP with lb and ub, both can be row/col vector
% doorvanbei
% 20220711
lbtest = lb>0;
ubtest = ub<inf;
lbnum = sum(lbtest);
ubnum = sum(ubtest);
lbpos = find(lbtest);
ubpos = find(ubtest);
l = length(c);
Allb = zeros(lbnum,l);
bllb = zeros(lbnum,1);
for i = 1:lbnum
    Allb(i,lbpos(i)) = -1;
    bllb(i) = -lb(lbpos(i));
end
Alub = zeros(ubnum,l);
blub = zeros(ubnum,1);
for i = 1:ubnum
    Alub(i,ubpos(i)) = 1;
    blub(i) = ub(ubpos(i));
end
[x,v,f] = myip(c,iind,[Al;Allb;Alub],[bl;bllb;blub],A,b);
function f = iscutok(lhs,rhs,ub)
r = abs(lhs);
f = 1;
lb = 1/ub;
if any( (r>0) & (r<lb) ) || any(r>ub)
    f = 0;
    return
end
r = abs(rhs);
ub = 10 * ub;
lb = 1/ub;
if ( (r>0) && (r<lb) ) || (r>ub)
    f = 0;
end
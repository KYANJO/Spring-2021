function c=minmod(a,b)
if(a*b<=0)
    c=0;
else
    if(abs(a)<abs(b))
        c=a;
    else
        c=b;
    end
end
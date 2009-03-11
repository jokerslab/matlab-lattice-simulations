function r=properLabel(Sn,LL)
    r=Sn;
    LL=[0;LL];
    t=r+1;
    t=-LL(t);
    while any(any(t>0))
        r(t>0)=t(t>0);
        q=-LL(t(t>0)+1);
        t(t>0)=q;
        
    end
end

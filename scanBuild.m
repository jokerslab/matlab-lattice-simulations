function list=scanBuild(list,pt,neighbors,s)
list=[list;pt];
n=neighbors(:,pt);
n=n(s(n)==s(pt));
if ~isempty(n)
    
    for i=1:length(n)
        if ~any(n(i)==list)
           list=scanBuild(list,n(i),neighbors,s);
        end
    end     
end
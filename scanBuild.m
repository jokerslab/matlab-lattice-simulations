function list=scanBuild(list,pt,neighbors,s)
list=[list;pt];
n=neighbors(:,pt);
%n=n(s(n)==s(pt));
n=n(s(n)>1);
while ~isempty(n)
    
    c=setdiff(n,list);
    list=[list;c];
    n=unique(neighbors(:,c));
    %n=n(s(n)==s(pt));
    n=n(s(n)>1);
end
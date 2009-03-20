function [c,LL,F]=mergeClusters(i,j,neighbor,c,LL,F)
%Takes unique cluster labels and merges them
[mn,imn]=min(neighbor);
                        
sumc = 1+sum(LL(neighbor));
sumx = j+sum(F.x(neighbor));
sumy = i+sum(F.y(neighbor));
F.x2(neighbor(imn))=j*j+sum(F.x2(neighbor));
F.y2(neighbor(imn))=i*i+sum(F.y2(neighbor));

%Label of labels negative sign points to
%'proper label'
LL(neighbor(neighbor~=mn))=-mn;



LL(neighbor(imn))=sumc;
F.x(neighbor(imn))=sumx;
F.y(neighbor(imn))=sumy;
c(i,j)=neighbor(imn);            
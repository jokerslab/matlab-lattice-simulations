clear;
s=ones(20);
s=double(rand(20)>.7)+1;
rows=20;
cols=20;
% s(1:5,1:5)=2;
% s(1:5,10:15)=2;
% s(10:15,1:5)=2;
% s(10:15,10:15)=2;
ne=zeros(8,20,20);
    for i=1:20
        for j=1:20
            ne(:,i,j)=neighbors2d(i,j,rows,cols);
        end
    end
[c labs csize cpos]=clusterCount(s,2);
plotLattice(c,labs)
bigcI=find(csize>1);
for t=1:100
 i=floor(rand*length(bigcI))+1;
    pts=find(c==labs(bigcI(i)));
    dir=floor(rand*4)+1;
    s(pts)=1;
    s(ne(dir,pts))=2;
    newpts=ne(dir,pts);
    newne=setxor(newpts,unique(ne(1:4,newpts)));
    

[c labs csize cpos]=clusterCount(s,2);
bigcI=find(csize>1);
    
    plotLattice(c, labs)
    pause(.2)
end
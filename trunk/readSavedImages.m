d=dir('./images');
e=dir(['./images/' d(35).name '/*.png']);

pts=length(e);
rel=zeros(pts,1);
relt=zeros(pts,1);
for i=1:length(e)
s=imread(['./images/' d(35).name '/' e(i).name])+1;
rel(i)=nnz(s==5);
relt(i)=500*(i-1);
end
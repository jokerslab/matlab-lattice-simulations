imagePerC=5;
d=dir('./images');
dc=[39];

for dt=1:length(dc)
e=dir(['./images/' d(dc(dt)).name '/*.png']);

pts=length(e);
rel{dt}=zeros(pts,1);

relt{dt}=zeros(pts,1);
larC{dt}=zeros(floor(pts/imagePerC)+1,1);
nnc{dt}=zeros(floor(pts/imagePerC)+1,1);
fre1{dt}=zeros(floor(pts/imagePerC)+1,1);
t=0;
for i=1:length(e)
    if mod(i-1,imagePerC)==0
        t=t+1;
        if mod(t,200)==0
            t
        end
        s=imread(['./images/' d(dc(dt)).name '/' e(i).name])+1;
        rel{dt}(i)=nnz(s==2);
        relt{dt}(i)=500*(i-1);
        [c LL f pts2]=clusterCountEHK2(s,[2 3]);
        larC{dt}(t)=max(LL);
        nnc{dt}(t)=nnz(LL>5);
        fre1{dt}(t)=nnz(LL==1);
    end
end
end
m=length(LL3);
mindex=1;
clusterBubbleX=zeros(m*20,1);
clusterBubbleY=zeros(m*20,1);
clusterBubbleS=zeros(m*20,1);
for i=1:m
    if (nnz(LL3{i}>0)>0)
        clusterSizes=unique(LL3{i}(LL3{i}>0));
        for j=1:length(clusterSizes)
            clusterBubbleX(mindex)=i;
            clusterBubbleY(mindex)=clusterSizes(j);
            clusterBubbleS(mindex)=10*nnz(LL3{i}(LL3{i}==clusterSizes(j)));
            mindex=mindex+1;
        end
        
    end
end
clusterBubbleX=clusterBubbleX(1:mindex);
clusterBubbleY=clusterBubbleY(1:mindex);
clusterBubbleS=clusterBubbleS(1:mindex);


crystalclusters=zeros(m,20);
for i=1:m
    for j=1:20
        crystalclusters(i,j)=nnz(LL3{i}>j);
            
        
    end
end
plot(y,fre1/4096,'b',1:1000000,frac3/4096,'r',1:1000000,frac2/4096,'g')
xlabel('Simulation Steps');
ylabel('Concentration')
legend('Free B','B*','B')
title('k_bT = 1 \chi_{AB} = -0.1 \chi_{AB^*} = 0.2 \epsilon = -0.5 \gamma= 0.5')
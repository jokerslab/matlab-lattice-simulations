%2d Lattice
%Ising Model
function [s,Energy,cor,clust,cave]=Ising2D(kbT,phi,s)

%rand('twister',sum(100*clock))
%randn('state',sum(100*clock))
%rand('twister',10)
%randn('state',5)
startTime=datestr(now);
Steps=20000;
Energy=zeros(Steps+1,1);
showplot=false;
saveplot=true;
Sps=1000;  %steps per save of data


if nargin<2,phi=.5; end
if nargin<1,kbT=1; end
disp('Setting up Initial state')
%s=sign(randn(rows,cols));
if nargin<3
    rows=50;
    cols=50;
    s=sign(randn(rows,cols)-sqrt(2)*erfcinv(2*phi));
else
    [rows,cols]=size(s);
end

N=rows*cols;

%4 neighbors in 2d
neighbors=zeros(4,rows,cols);

%how many 1's do we have
m=sum(sum(s==1));
plotLattice(s)
title('Initial State')
if saveplot
filestr=['nucMC/pics/' startTime 'Initialpos_T_' num2str(kbT) 'phi_' num2str(phi)];
print('-r600','-depsc',filestr);
close
end
%get the possible energy types
d_en=testState([-1 1 2]); 

%Classify pair types;
[px Energy(1)]=classifypairs2(s);

%Create Neighbor lists
for i=1:rows
    for j=1:cols
        neighbors(:,i,j)=neighbors2(i,j,rows,cols);
    end
end



%Create lists for each scenerio pair
 
for i=1:length(d_en)
    p{i}=find(px==d_en(i));
end
%Create energy for each sceneriio pair
%ignore first value in d_en since it indicates swapping same state
d_energy=[0; d_en(2:end)];


d_p=exp(-d_energy./kbT);

cor(:,:,1)=LatticeCorrelation(s,5);
s_old=s;
[c,labels,csize]=clusterCount(s);
clust(1)=sum(csize>10);
cave(1)=mean(csize(csize>10));
disp('Starting MainLoop')
%Number of steps 
oldavE=Energy(1);
for t=1:Steps
    if mod(t,Sps)==0
        t
        %LatticeCorrelation(s,0,s_old)
        avE=mean(Energy(t-Sps+1:t));
        avE-oldavE
        oldavE=avE;
        s_old=s;
        cor(:,:,t/Sps+1)=LatticeCorrelation(s,5);
        [c,labels,csize]=clusterCount(s);
        clust(t/Sps+1)=sum(csize>10);
        cave(t/Sps+1)=mean(csize(csize>10));
        
    end
        
    %Running sum of number of pairs * delta energy of switch
    %Skip the first pair category since the swap does not change
    q1(1)=length(p{2})*d_p(2);
    for r=2:length(d_p)-1
        q1(r)=q1(r-1)+length(p{r+1})*d_p(r+1);
    end
    
    %Select a random number to determine which category pair to swap
       temp=rand*q1(end);
       st=find(q1>temp,1)+1;
       %Select a random pair from that category
       q=p{st}(ceil(length(p{st}).*rand));
       Energy(t+1)=Energy(t)+d_energy(st);
       %Now that q=the pair index lets do the swap and recalculate the
       %list
    if mod(q,2)==0  %Horizontal bond
        si=q/2;   %bonds 1 & 2 belong to site 1 ...
        nd=4;     %neighbor 4 is to the 'right'
    else  %Vertical bond
        si=(q+1)/2;
        nd=2;    %neighbor 2 is down
    end
        s(si)=-s(si);  %"swap" values 2 state model can just reverse sign
        
        %get the row and column of the site
        [r c]=ind2sub([rows cols],si);
        if showplot
        line('xdata',c,'ydata',-r,'markersize',10,'marker','o', ...
    'linestyle','none','markerfacecolor','red','markeredgecolor','black')
        end
        %find the neighbors of the site
        n_1=neighbors(:,r,c);
        
        %Swap the neighbor of the site that corresponds to the pair
        s(n_1(nd))=-s(n_1(nd));
        %for now, reclassify bonds
        %in the future can precalculate changes
        
        %find the rows and columns of the neighbors of other pair
        [r c]=ind2sub([rows cols],n_1(nd));
        if showplot
        line('xdata',c,'ydata',-r,'markersize',10,'marker','o', ...
    'linestyle','none','markerfacecolor','red','markeredgecolor','black')
        end
        n_1b=[n_1; neighbors(:,r,c)];
        [r c]=ind2sub([rows cols],n_1b);
        
        %find the neighbors of the neighbors ...
        %these pair energies are 
        %afftect by the swap
        n_2=[];
        for jj=1:length(r)
            n_2=[n_2; neighbors(:,r(jj),c(jj))];
        end
        nn=unique([n_1b; n_2]);
        [r1 c1]=ind2sub([rows cols],nn);
        
        
        %Loop through the pairs from all the neighbors and 
        %reclassify the type
        %update the lists if things change
        [p2,e0,pI]=classifypairs2(s,r1,c1);
        for kk=1:length(pI)
            if px(pI(kk))~=p2(pI(kk))
                oldgroup=find(d_en==px(pI(kk)));
                newgroup=find(d_en==p2(pI(kk)));
                p{oldgroup}(p{oldgroup}==pI(kk))=[];
                p{newgroup}=[p{newgroup}; pI(kk)];
                
            end
        end
        px(pI)=p2(pI);
    
%     for kk=1:length(r1)
%             i=r1(kk);
%             ip1=i+1;
%             if i==rows
%                 ip1=1;
%             end
%             j=c1(kk);
%             jp1=j+1;
%               if j==cols
%                 jp1=1;
%               end
%               oldp=pairs(2*i-1,j);
%             if s(i,j)==s(ip1,j)
%                 pairs(2*i-1,j)=-1;
%             else
%                 s1=sum(s(i,j)==s(neighbors2(i,j,rows,cols)));
%                 s2=sum(s(ip1,j)==s(neighbors2(ip1,j,rows,cols)));
%                 pairs(2*i-1,j)=s1+s2;
%             
%             end
%             if oldp~=pairs(2*i-1,j)
%                 ind=(j-1)*2*rows+2*i-1;
%                 p{oldp+2}((p{oldp+2}==ind))=[];
%                 p{pairs(2*i-1,j)+2}=[p{pairs(2*i-1,j)+2};ind];
%                 
%             end
%             oldp=pairs(2*i,j);
%             if s(i,j)==s(i,jp1)
%                 pairs(2*i,j)=-1;
%             else
%                 s1=sum(s(i,j)==s(neighbors2(i,j,rows,cols)));
%                 s2=sum(s(i,jp1)==s(neighbors2(i,jp1,rows,cols)));
%                 pairs(2*i,j)=s1+s2;
%             end
%             if oldp~=pairs(2*i,j)
%                 ind=(j-1)*2*rows+2*i;
%                 p{oldp+2}((p{oldp+2}==ind))=[];
%                 p{pairs(2*i,j)+2}=[p{pairs(2*i,j)+2};ind];
%                 
%             end
%             
%     end
       if showplot
       %pause
       subplot(3,3,[1:6])
       plotLattice(s)
       plotBonds(px,d_en(2:end))
       pause(.01)
       subplot(3,3,7)
       for i=2:length(p)
           ng(i-1)=length(p{i});
       end
       bar(d_en(2:end),ng)
       subplot(3,3,8)
       plot(0:t,Energy(1:t+1));
       subplot(3,3,[1:6])
       %pause
       end
end
figure()
plotLattice(s);
title(['State after ' num2str(t) ' steps'])
if saveplot
    filestr=['nucMC/pics/' startTime 'Finalpos'];
    print('-r600','-depsc',filestr);
end
%figure()
%plotLattice(pairs,[0 1 2 3 4 5 6 -1])
figure()
%plotLattice(classifypairs(s),[0 1 2 3 4 5 6 -1])
plot(0:t,Energy(1:t+1));
ylabel('Energy')
xlabel('Steps')
title('Energy of System')
if saveplot
filestr=['nucMC/pics/' startTime 'Energy'];
print('-r600','-depsc',filestr);
end
ts=0:Sps:t;
%figure
%plot(ts,cor(5,6,:))
%title('Correlation 1 lattice site away')
%filestr=['nucMC/pics/' startTime 'Correlation'];
%print('-r600','-depsc',filestr);
figure
plot(ts,cave)
title('Average cluster (>10) size')
if saveplot
filestr=['nucMC/pics/' startTime 'Csize'];
print('-r600','-depsc',filestr);
end
end






function pairs=classifypairs(s)
%Classify the bond pair by adding up the number of like-like bonds 
%among the neighbors
[rows,cols]=size(s);
pairs=zeros(2*rows,cols)+25;
for i=1:rows
    ip1=i+1;
    if i==rows
        ip1=1;
    end
    for j=1:cols
       jp1=j+1;
       if j==cols
           jp1=1;
       end
        if s(i,j)==s(ip1,j)
            pairs(2*i-1,j)=-1;
        else
            s1=sum(s(i,j)==s(neighbors2(i,j,rows,cols)));
            s2=sum(s(ip1,j)==s(neighbors2(ip1,j,rows,cols)));
            pairs(2*i-1,j)=s1+s2;
            
        end

        if s(i,j)==s(i,jp1)
            pairs(2*i,j)=-1;
        else
            s1=sum(s(i,j)==s(neighbors2(i,j,rows,cols)));
            s2=sum(s(i,jp1)==s(neighbors2(i,jp1,rows,cols)));
            pairs(2*i,j)=s1+s2;
        end
    end
end
end

function n=neighbors2(i,j,rows,columns)
        index=(j-1)*rows+i;
        n=zeros(1,4);
        if (i~=1) && (i~=rows)
            n(1:2)=[index-1 index+1];
        elseif i==1
            n(1:2)=[index+rows-1 index+1];
        else
            n(1:2)=[index-1 index-rows+1];
        end
        if (j~=1) && (j~=columns)
            n(3:4)=[index-rows index+rows];
        elseif j==1
            n(3:4)=[index+(columns-1)*rows index+rows];
        else
            n(3:4)=[index-rows index-(columns-1)*rows];
        end
        
    end
function n=pbc(i,L)
    n=i;
    if i<1 
        n=i+L;
    end
    if i>L
        n=i-L;
    end
end
%Ising Model
function [s,cs,LL,F,Energy]=Ising2D(kbT,phi,s,Steps)

%Random Number Generator Setup
%rand('twister',sum(100*clock))
%randn('state',sum(100*clock))
%rand('twister',10)
%randn('state',5)


%Bond energies
SS=1;
SA=2;
SC=-2;
AA=1;
AC=-1;
CC=4;
et=[SS SA SC;SA AA AC;SC AC CC];

NUMBER_NEIGHBORS=8;
states=[1 2 3];
%1--solvent
%2 --amorphous
%3 --Crystal

%Set up Default Arguments

if nargin<2,phi=.5; end
if nargin<1,kbT=1; end
disp('Setting up Initial state')
%s=sign(randn(rows,cols));
if nargin<3
    %Standard size 25,25
    rows=25;
    cols=25;
    s=sign(randn(rows,cols)-sqrt(2)*erfcinv(2*phi));
    s(s==1)=2;
    s(s==-1)=1;
else
    %If we input the size of the grid
    if length(s)==2
        rows=s(1);
        cols=s(2);
        s=sign(randn(rows,cols)-sqrt(2)*erfcinv(2*phi));
        s(s==1)=2;
        s(s==-1)=1;
    else
    [rows,cols]=size(s);
    end
end
if nargin<4
    Steps=1000;
end

%Output Parameters
Sps=1000;  %steps per save of data
showplot=true;
saveplot=false;
saveFrames=false;
startTime=datestr(now);
animate=false;
filestr=['home/pjung/nucMC/pics/' startTime '_T_' num2str(kbT) 'phi_' num2str(phi) '_'];




N=rows*cols;

%Initialize Arrays
Energy=zeros(Steps+1,1);

%4 neighbors in 2d
%8 next nearest neighbors in 2d
neighbors=zeros(NUMBER_NEIGHBORS,rows,cols);

%how many 1's do we have
m=sum(sum(s==1));


plotLattice(s,states)
title('Initial State')

if saveplot
print('-r600','-depsc',[filestr 'InitialPos']);
close
end


%get the possible energy types
%d_en=testState([-1 1 2]); 


%Create Neighbor lists
for i=1:rows
    for j=1:cols
        neighbors(:,i,j)=neighbors2d(i,j,rows,cols);
    end
end
%Index of Neighbors
%xxx516xxx
%xxx3o4xxx
%xxx728xxx

%Not yet fully implemented
%need to work out details to classify possible moves
swapMoves=[2 4 5 6];
flipMoves=[2];
nMoves=length(swapMoves)+length(flipMoves);

%Classify pair types;
[px Energy(1)]=classifypairs2(s,et,swapMoves,neighbors);
d_en=unique(px);
%Create lists for each scenerio pair
 
for i=1:length(d_en)
    p{i}=find(px==d_en(i));
end

plistindex=zeros(size(px));
for i=1:length(p)
    for j=1:length(p{i})
        plistindex(p{i}(j))=j;
    end
end
%Create energy for each sceneriio pair
%ignore first value in d_en since it indicates swapping same state
d_energy=[0; d_en(2:end)];
zerogroup=find(d_en==0);
if isempty(zerogroup)
        zerogroup=length(p)+1;
        d_en(zerogroup)=0;
        d_energy(zerogroup)=0;
        p{zerogroup}=[];
        
end
d_p=exp(-d_energy./kbT);

cor(:,:,1)=LatticeCorrelation(s,5);
%s_old=s;
[cs,LL,F]=clusterCountEHK2(s,2);
cs=properLabel(cs,LL);
largest_label=length(LL);
%maxLabel=max(labels)+1;
%clust(1)=sum(csize>10);
%cave(1)=mean(csize(csize>10));

disp('Starting MainLoop')
%Number of steps 
oldavE=Energy(1);
for t=1:Steps
    
    
    %Intermittently save the current state for output
    if mod(t,Sps)==0
        t
        %LatticeCorrelation(s,0,s_old)
        %avE=mean(Energy(t-Sps+1:t));
        %avE-oldavE
        %oldavE=avE;
        
        cor(:,:,t/Sps+1)=LatticeCorrelation(s,5);
        %[c,labels,csize]=clusterCount(s,3);
        %clust(t/Sps+1)=sum(csize>10);
        %cave(t/Sps+1)=mean(csize(csize>10));
        
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
%        si=mod(q,N);
%        if si==0,si=N;end
%        nd=nmoves(floor((q-1)/N)+1);
if q<0
    pts=find(properLabel(cs,LL)==-q);
    dir=floor(rand*4)+1;
    s(pts)=1;
    s(neighbors(dir,pts))=2;
    newpts=neighbors(dir,pts);
    pts=unique([pts;newpts']);
    nn=unique(neighbors(:,neighbors(:,pts)));
    [cs,LL,F]=clusterCountEHK2(s,[2 3]);
    cs=properLabel(cs,LL);
    largest_label=length(LL);
else
       [nd,si]=ind2sub(size(px),q);
       if nd<=length(swapMoves)
           
        nd=swapMoves(nd);
        temp=s(si);
       
        n_1=neighbors(:,si);
        s(si)=s(n_1(nd));  %"swap" values
        
        %Swap the neighbor of the site that corresponds to the pair
        s(n_1(nd))=temp;
        
        updateClusters
        %get the row and column of the site
       else
           %we flip states
           nd=1;
           if s(si)==2
                s(si)=3;
           else
               s(si)=2;
               
           end
           n_1=si;
           
       end
        if showplot
            [r c]=ind2sub([rows cols],si);
        line('xdata',c,'ydata',-r,'markersize',10,'marker','o', ...
    'linestyle','none','markerfacecolor','red','markeredgecolor','black')

            [r c]=ind2sub([rows cols],n_1(nd));
        line('xdata',c,'ydata',-r,'markersize',10,'marker','o', ...
    'linestyle','none','markerfacecolor','red','markeredgecolor','black')
        
        end
        
        %for now, reclassify bonds
        %in the future can precalculate changes
        
        %find the rows and columns of the neighbors of other pair
        
        
        n_1b=[n_1; neighbors(:,n_1(nd))];
        
        %find the neighbors of the neighbors ...
        %these pair energies are 
        %afftect by the swap
       
        nn=unique(neighbors(:,n_1b));
        
        
        %Loop through the pairs from all the neighbors and 
        %reclassify the type
        %update the lists if things change
end
        [p2,e0]=classifypairs2(s,et,swapMoves,neighbors,nn);
        
         p{zerogroup}(p{zerogroup}<0)=[];
        
        for kk=1:length(nn)
            for pp=1:nMoves

                if px(pp,nn(kk))~=p2(pp,nn(kk))
                    oldgroup=find(d_en==px(pp,nn(kk)));
                    newgroup=find(d_en==p2(pp,nn(kk)));
                    %The following search is slow
                    ind=plistindex(pp,nn(kk));
                    
                    
                        plistindex(p{oldgroup}(ind:end))=...
                        plistindex(p{oldgroup}(ind:end))-1;
                    
                    p{oldgroup}(ind)=[];
                    if isempty(newgroup)
                        newgroup=length(p)+1;
                        d_en(newgroup)=p2(pp,nn(kk));
                        d_energy(newgroup)=p2(pp,nn(kk));
                        d_p=exp(-d_energy./kbT);
                        p{newgroup}=[nMoves*(nn(kk)-1)+pp ];
                    else
                        
                        p{newgroup}=[p{newgroup};nMoves*(nn(kk)-1)+pp ];
                        
                    end
                    
                    
                    
                    plistindex(pp,nn(kk))=length(p{newgroup});
                    
                end
            end
        end
        px(:,nn)=p2(:,nn);
%         for kk=1:length(pI)
%             if px(pI(kk))~=p2(pI(kk))
%                 oldgroup=find(d_en==px(pI(kk)));
%                 newgroup=find(d_en==p2(pI(kk)));
%                 %The following search is slow
%                 %should keep track of the index
%                 p{oldgroup}(p{oldgroup}==pI(kk))=[];
%                 p{newgroup}=[p{newgroup}; pI(kk)];
%                 
%             end
%         end
%         px(pI)=p2(pI);
        
%move clusters
        %[c labs csize cpos]=clusterCount(s,2);
        bigcI=find(LL>1);
   
%    p{zerogroup}=[p{zerogroup};-bigcI];
    
        
        
       if showplot
       if strcmp(get(gcf,'CurrentCharacter'),'p')    
        pause
       end
       
       %subplot(3,3,[1:6])
       plotLattice(s, states)
       %plotBonds(px,d_en(2:end))
       title(['Step' num2str(t,'%05d')]);
       pause(.01)
       %subplot(3,3,7)
       %for i=2:length(p)
       %    ng(i-1)=length(p{i});
       %end
       %bar(d_en(2:end),ng)
       %subplot(3,3,8)
       %plot(0:t,Energy(1:t+1));
       %subplot(3,3,[1:6])
       if saveFrames
%         mF=getframe(gcf);
%         mFI=frame2im(mF);
         directory = 'images/';  % The next four lines parse and assemble fil
%     
         filename = [directory, num2str(t,'%06d'),'.png'];
% 
%         imwrite(mFI,filename, 'png'); % Finally write individual images
%         clear mF
%         clear mFI
            set(gcf,'PaperPositionMode','auto')
            print(gcf,'-dpng','-r0',filename)
        end
       if get(gcf,'CurrentCharacter')=='p'    
        pause
       end
       end
end

%Plot the outputs

figure()
plotLattice(s, [1 2 3],{'none' [.5 .5 .5],'black'})
title(['State after ' num2str(t) ' steps'])
if saveplot
    
    print('-r600','-depsc',[filestr 'Finalpos']);
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
%figure
%plot(ts,cave)
%title('Average cluster (>10) size')
%if saveplot
%filestr=['nucMC/pics/' startTime 'Csize'];
%print('-r600','-depsc',filestr);
%end
end
% function pausefig(src,evnt)
%       if evnt.Character == 't'
%          pause(1)
%       end
% end


function [s,cs,LL,F,Energy]=LatticeSimulation(varargin)
%Lattice Simulation Model main routine
%This function is the main method used to run simulations
%This is the main function/routine used to start and compute the simulation. 
%There are several available parameters and options that can be set,
%although the simulation will run with 'default' arguments if no parameters
%are set
%
%List of Parameters
%
%kbT    Temperature in units of kbT (default is 1)
%phi    Volume fraction of solute if generating random initial state (0.5)
%Steps  # of simulation steps to run (default is 1000)
%
%List of options
%Examples
%Running default simulation:
%LatticeSimulation();
%
%Running default simulation but with 1 million step simulation and kbt=2:
%LatticeSimulation('kbt',2,'Steps',1000000);

%LatticeSimulation(kbT,phi,s,Steps,chiAB,chiAC,epsilon,gamma)

% Create an instance of the inputParser class.
p = inputParser;
p.addParamValue('kbT',1,@(x)isnumeric(x) && isscalar(x) && x>0);
p.addParamValue('phi',.5,@(x)isnumeric(x) && isscalar(x) && x>=0 &&x<=1);
p.addParamValue('s',[64 64]);
p.addParamValue('Steps',1000,@(x)isnumeric(x) && isscalar(x));
p.addParamValue('chiAB',-0.1,@(x)isnumeric(x) && isscalar(x));
p.addParamValue('chiAC',0.6,@(x)isnumeric(x) && isscalar(x));
p.addParamValue('epsilon',-1,@(x)isnumeric(x) && isscalar(x));
p.addParamValue('gamma',1.1,@(x)isnumeric(x) && isscalar(x));

%Output Parameters
p.addParamValue('Sps',100,@(x)isnumeric(x) && isscalar(x)); %steps per save of data
p.addParamValue('Spp',3000,@(x)isnumeric(x) && isscalar(x)); %steps per plot
p.addParamValue('Spf',1500,@(x)isnumeric(x) && isscalar(x)); %steps per frame used in creating movies
p.addParamValue('showplot',true,@(x)islogical(x) && isscaler(x));
p.addParamValue('saveplot',false,@(x)islogical(x) && isscaler(x));
p.addParamValue('saveFrames',true,@(x)islogical(x) && isscaler(x));
p.addParamValue('showMoves',false,@(x)islogical(x) && isscaler(x));

p.parse(varargin{:});
% Show the value of a specific argument.
params=p.Results;
clear p;

%Random Number Generator Setup
%rand('twister',sum(100*clock))
%randn('state',sum(100*clock))
%rand('twister',10)
%randn('state',5)

    
et=computeInteractionParameters(params.chiAB,params.chiAC,params.epsilon,params.gamma);

NUMBER_NEIGHBORS=8;

states=[1 2 3 4 5];
%1--solvent
%2 --amorphous
%3 --Crystal
%4 --Wall
%5 --Water

%Set up Default Arguments


disp('Setting up Initial state')


    %If we input the size of the grid
    if length(params.s)==2
        rows=params.s(1);
        cols=params.s(2);
        s=createInitialState(params.phi,[rows cols]);
    else
        [rows,cols]=size(params.s);
        s=params.s;
    end



startTime=datestr(now,'yyyymmmdd-HHMM');
animate=false;

if params.saveFrames
    mkdir(['images/' startTime]);
    directory = strcat('images/',startTime,'/'); 
    filename = [directory 'Em_T_' num2str(params.kbT) '.txt'];
    save(filename,'et','-ASCII')
    filename = [directory, num2str(0,'%07d'),'.png'];
    imwrite(s,[1 1 1;0 0 0;0 1 0;.5 .5 .5; 0 0 1],filename)
end

orientation=zeros(size(s));

checkAbsorb=all(s(1,:)==5);

N=rows*cols;

%Initialize Arrays
Energy=zeros(params.Steps+1,1);
frac2=zeros(params.Steps,1);
frac3=zeros(params.Steps,1);

%4 neighbors in 2d
%8 next nearest neighbors in 2d
neighbors=zeros(NUMBER_NEIGHBORS,rows,cols);

%how many 1's do we have
m=sum(sum(s==1));


plotLattice(s,states)
title('Initial State')

if params.saveplot
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
global pairsMap
pairsMap=-1000*ones(3,3,3,3,3,3,3,3,3,3,3,3);
[px Energy(1)]=classifypairs2(s,et,swapMoves,neighbors);
d_en=unique(round(px*1e5)/1e5);
%Create lists for each scenerio pair
 
for i=1:length(d_en)
    p{i}=(abs(px-d_en(i))<1e-5);
    pLength(i)=nnz(p{i});
end


%Create energy for each sceneriio pair
%ignore first value in d_en since it indicates swapping same state
d_energy=[0; d_en(2:end)];
zerogroup=find(d_en==0);
if isempty(zerogroup)
        zerogroup=length(p)+1;
        d_en(zerogroup)=0;
        d_energy(zerogroup)=0;
        p{zerogroup}=false(size(px));
        pLength(zerogroup)=0;
        
end
d_p=exp(-d_energy./params.kbT);

cor(:,:,1)=LatticeCorrelation(s,5);
%s_old=s;
[cs,LL,F]=clusterCountEHK2(s,[2 3]);
cs=properLabel(cs,LL);
largest_label=length(LL);
larC(1)=max(LL);
nnc(1)=nnz(LL>5);
fre1(1)=nnz(LL==1);
%maxLabel=max(labels)+1;
%clust(1)=sum(csize>10);
%cave(1)=mean(csize(csize>10));

disp('Starting MainLoop')
%Number of steps 
oldavE=Energy(1);
release=0;
sttt=1;
relt=0;

for t=1:params.Steps
    
    %Absorbing top boundry
    if checkAbsorb
    atb=(s(1,:)==2);
    if any(atb)
        sttt=sttt+1;
        release(sttt)=release(sttt-1)+nnz(atb);
        relt(sttt)=t;
        s(1,atb)=5;
    end
    end
    %Intermittently save the current state for output
    frac2(t)=nnz(s==2);
    frac3(t)=nnz(s==3);
    if mod(t,params.Sps)==0
        t
        %LatticeCorrelation(s,0,s_old)
        %avE=mean(Energy(t-Sps+1:t));
        %avE-oldavE
        %oldavE=avE;
        [cs,LL,F]=clusterCountEHK2(s,[2 3]);
        LL23{t/params.Sps+1}=LL;
        larC(t/params.Sps+1)=max(LL);
        nnc(t/params.Sps+1)=nnz(LL>5);
        fre1(t/params.Sps+1)=nnz(LL==1);
        [cs,LL,F]=clusterCountEHK2(s,[3]);
        LL3{t/params.Sps+1}=LL;
        %cor(:,:,t/Sps+1)=LatticeCorrelation(s,5);
        %[c,labels,csize]=clusterCount(s,3);
        %clust(t/Sps+1)=sum(csize>10);
        %cave(t/Sps+1)=mean(csize(csize>10));
        
    end
        
    %Running sum of number of pairs * delta energy of switch
    %Skip the first pair category since the swap does not change
    q1(1)=pLength(2)*d_p(2);
    for r=2:length(d_p)-1
        q1(r)=q1(r-1)+pLength(r+1)*d_p(r+1);
    end
    
    %Select a random number to determine which category pair to swap
       if q1(end)==0
           break;
       end
       temp=rand*q1(end);
       st=find(q1>temp,1)+1;
       %Select a random pair from that category
       tList=find(p{st});
       q=tList(ceil(pLength(st).*rand));
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
       n_1=neighbors(:,si);
       if nd<=length(swapMoves)
           
        nd=swapMoves(nd);
        temp=s(si);
       
        
        if temp==2 && s(n_1(nd))==2 %flip bond
            s(si)=3;
            s(n_1(nd))=3;
            orientation(si)=nd;
            orientation(n_1(nd))=nd;
        else
          s(si)=s(n_1(nd));  %"swap" values
        
          %Swap the neighbor of the site that corresponds to the pair
          s(n_1(nd))=temp;
        end
%        updateClusters
        %get the row and column of the site
       else
           %we flip states
           
           nd=1;
           if s(si)==2
                s(si)=3;
           else
               s(si)=2;
                for gi=1:4
                    if s(neighbors(gi,si))==3
                        
                        if ~any(s(neighbors(1:4,neighbors(gi,si)))==3)
                            s(neighbors(gi,si))=2;
                        end
                    end
                  
                end
               
           end
          % n_1=si;
           
       end
       
       
       %No longer need to look at moves indivdually
       
       
%         if showplot && showMoves
%             [r c]=ind2sub([rows cols],si);
%         line('xdata',c,'ydata',-r,'markersize',10,'marker','o', ...
%     'linestyle','none','markerfacecolor','red','markeredgecolor','black')
% 
%             [r c]=ind2sub([rows cols],n_1(nd));
%         line('xdata',c,'ydata',-r,'markersize',10,'marker','o', ...
%     'linestyle','none','markerfacecolor','red','markeredgecolor','black')
%         
%         end
        
        


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
        
%         [crap1 EnergyTemp]=classifypairs2(s,et,swapMoves,neighbors);
%         if abs(EnergyTemp-Energy(t+1))>.01
%            EnergyTemp
%            Energy(t+1)
%         end
        % p{zerogroup}(p{zerogroup}<0)=[];
        
        for kk=1:length(nn)
            for pp=1:nMoves

                if abs(px(pp,nn(kk))-p2(pp,nn(kk)))>0.001
                    oldgroup=find(abs(d_en-px(pp,nn(kk)))<0.001);
                    newgroup=find(abs(d_en-p2(pp,nn(kk)))<0.001);
                    
                    %Try using logical trackers
                    p{oldgroup}(pp,nn(kk))=0;
                    pLength(oldgroup)=pLength(oldgroup)-1;
                    %The following search is slow
                    
%                     ind=plistindex(pp,nn(kk));
%                     
%                     
%                         plistindex(p{oldgroup}(ind:end))=...
%                         plistindex(p{oldgroup}(ind:end))-1;
%                     
%                     p{oldgroup}(ind)=[];
                    if isempty(newgroup)
                        newgroup=length(p)+1;
                        d_en(newgroup)=p2(pp,nn(kk));
                        d_energy(newgroup)=p2(pp,nn(kk));
                        d_p=exp(-d_energy./params.kbT);
                        %p{newgroup}=[nMoves*(nn(kk)-1)+pp ];
                        p{newgroup}=false(size(px));
                        pLength(newgroup)=0;
                        
                    else
                        
                        %p{newgroup}=[p{newgroup};nMoves*(nn(kk)-1)+pp ];
                        
                    end
                    
                        p{newgroup}(pp,nn(kk))=1;
                        pLength(newgroup)=pLength(newgroup)+1;
                    
                    %plistindex(pp,nn(kk))=length(p{newgroup});
                    
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
        %bigcI=find(LL>1);
   
%    p{zerogroup}=[p{zerogroup};-bigcI];
    
        
        
       if params.showplot  && mod(t,params.Spp)==0
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
       
       if get(gcf,'CurrentCharacter')=='p'    
        pause
       end
       end
       if params.saveFrames && mod(t,params.Spf)==0

          % The next four lines parse and assemble fil
%     
         filename = [directory, num2str(t,'%07d'),'.png'];
% 
%            set(gcf,'PaperPositionMode','auto')
            
%            print(gcf,'-dpng','-r0',filename)
            
            imwrite(s,[1 1 1;0 0 0;0 1 0;.5 .5 .5;0 0 1],filename)
            
        end
end

%Plot the outputs
Ising2Doutput



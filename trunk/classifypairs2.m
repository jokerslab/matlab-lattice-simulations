function [pairs,currentEnergy]=classifypairs2(state,energyMatrix,moves,ne,pts)
%Classify the bond pair by calculating change in energy if swapped

%pairs stores energy for swaping sites
%not swappable value to use
nomove=-250;


%flips stores energy difference for allowed flips

s=state;
% nStates=length(stateValues);
% logicalIndex=false([size(s) nStates]);
% for istate=1:nStates
%     logicalIndex(:,:,istate)=(s==stateValues(istate));
% end
% for istate=1:nStates
%     s(logicalIndex(:,:,istate))=istate;
% end
et=energyMatrix;
[rows,cols]=size(s);

if nargin<4
    ne=zeros(8,rows,cols);
    for i=1:rows
        for j=1:cols
            ne(:,i,j)=neighbors2d(i,j,rows,cols);
        end
    end
end
if nargin<5
    doFull=true;
    pts=1:rows*cols;
else
    doFull=false;
end


pairs=zeros(5,rows,cols);

currentEnergy=0;

for kk=1:length(pts)
    pts1=ne([2 4 5 6],pts(kk));
        if s(pts(kk))==1 || s(pts(kk))==4
            pairs(5,pts(kk))=nomove;
        end
        if s(pts(kk))==2
             pairs(5,pts(kk))=findEnergyDifference1(pts(kk),3);
        end
         if s(pts(kk))==3
               pairs(5,pts(kk))=findEnergyDifference1(pts(kk),2);
         end
    for cnt=1:length(pts1)
%         pairs(cnt+4,pts(kk))=nomove;
%         pairs(cnt+8,pts(kk))=nomove;

       
        
    if s(pts(kk))==s(pts1(cnt))  || s(pts(kk))==4 || s(pts1(cnt))==4
           pairs(cnt,pts(kk))=nomove;
    else
        pairs(cnt,pts(kk))=findEnergyDifference(pts(kk),pts1(cnt),s(pts1(cnt)),s(pts(kk)));
            
    end
    end
        

       
        if doFull
            currentEnergy=currentEnergy-et(s(pts(kk)),s(pts1(1)))-et(s(pts(kk)),s(pts1(2)));
        end
   
end
function Ediff=findEnergyDifference(p1,p2,newp1,newp2)
oldp1=s(p1);
oldp2=s(p2);
e11=sum(et(s(p1),s(ne(1:4,p1))));
e12=sum(et(s(p2),s(ne(1:4,p2))));
            %Swap
            s(p1)=newp1;
            s(p2)=newp2;
            %Energy After swap
e21=sum(et(s(p1),s(ne(1:4,p1))));
e22=sum(et(s(p2),s(ne(1:4,p2))));
            Ediff=-(e21+e22-e11-e12);
            %Swap Back for further calculations
            
s(p1)=oldp1;
s(p2)=oldp2;

end
function Ediff=findEnergyDifference1(p1,newp1)
oldp1=s(p1);

e11=sum(et(s(p1),s(ne(1:4,p1))));

            %Swap
            s(p1)=newp1;
           
            %Energy After swap
e21=sum(et(s(p1),s(ne(1:4,p1))));

            Ediff=-(e21-e11);
            %Swap Back for further calculations
            
s(p1)=oldp1;


end
end

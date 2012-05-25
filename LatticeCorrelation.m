function [g]=LatticeCorrelation(s,n,s2)
%Lattice Correlation 
%Calculate the autocorrelation of the lattice as it is shifted in space
[r,c]=size(s);
N=r*c;
c1=[];
if nargin<3
    s2=s;
end
g=zeros(2*n+1);
if n<=0
    g=sum(sum(s.*s2))/N;
else
for i=-n:n
    for j=-n:n
        if i>=0 && j>=0
            s3=[s2(i+1:r,j+1:c) s2(i+1:r,1:j);s2(1:i,j+1:c) s2(1:i,1:j)];
            
        end
        if i>=0 && j<0
            s3=[s2(i+1:r,c+j+1:c) s2(i+1:r,1:c+j);s2(1:i,c+j+1:c) s2(1:i,1:c+j)];
            
        end
        if i<0 && j>=0
            s3=[s2(r+i+1:r,j+1:c) s2(r+i+1:r,1:j);s2(1:r+i,j+1:c) s2(1:r+i,1:j)];
            
        end
        if i<0 && j<0
            s3=[s2(r+1+i:r,c+j+1:c) s2(r+1+i:r,1:c+j);s2(1:r+i,c+j+1:c) s2(1:r+i,1:c+j)];
            
        end
        g(i+n+1,j+n+1)=sum(sum(s.*s3))/N;
                
    end
end
end
%if n>0
%     su=[s2(2:r,:);s2(1,:)];
%     sd=[s2(r,:);s2(1:r-1,:)];
%     sl=[s2(:,2:c) s2(:,1)];
%     sr=[s2(:,c) s2(:,1:c-1)];
%     p1=sum(sum(abs(su-s)))/N;
%     p2=sum(sum(abs(sd-s)))/N;
%     p3=sum(sum(abs(sl-s)))/N;
%     p4=sum(sum(abs(sr-s)))/N;
%     c1=max([p1 p2 p3 p4]);
%     d1=LatticeCorrelation(s,n-1,su);
%     d2=LatticeCorrelation(s,n-1,sd);
%     d3=LatticeCorrelation(s,n-1,sl);
%     d4=LatticeCorrelation(s,n-1,sr);
%     c2=max([d1 d2 d3 d4]);
%     
%     c1=[c1 LatticeCorrelation(s,n-1,sr)];
    %g=[sum(sum(su.*s))/N;sum(sum(sd.*s))/N;sum(sum(sl.*s))/N;sum(sum(sr.*s))/N];
    %pick a mov
%end
    
    
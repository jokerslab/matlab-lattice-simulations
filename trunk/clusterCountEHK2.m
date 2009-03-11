function [c LL F]=clusterCountEHK2(s,value)
if nargin<2
    value=-1;
end
Not=value+25;
[r,c1]=size(s);
c=zeros(size(s));
largest_label=0;
LL=zeros(r*c1,1);
F=struct('x',cell(1,r*c1),'y',cell(1,r*c1));
neighbor=zeros(4,1);
for i=1:r
    im1=i-1;   
    for j=1:c1
        jm1=j-1;
        if s(i,j)==value;

            if i==1
                above=Not;
            else
                above=s(im1,j);
            end
            
            if j==1
                left=Not;
            else
                left=s(i,jm1);
            end
            if j==c1
                right=s(i,1);
            else
                right=Not;
            end
            if i==r
                below=s(1,j);
            else
                below=Not;
            end
            
                neighbor=zeros(4,1);
                
                if above==value
                    pl=properLabel(c(im1,j));
                   neighbor(1)=pl;
                end
                if left==value
                    neighbor(2)=properLabel(c(i,jm1));
                    
                end
                if right==value
                    neighbor(3)=properLabel(c(i,1));
                    
                end
                if below==value
                    neighbor(4)=properLabel(c(1,j));
                    
                end
               
                   if any(neighbor>0)
                        % NEED to UPDATE and adjust 
                        % x and y if it crosses lattice
                       [mn,imn]=min(neighbor(neighbor>0));
                        sum=1;
                        
                        sumx=j;
                        
                        sumy=i;
                        for ni=1:numneighbors
                            sum=sum+LL(neighbor(ni));
                            sumx=sumx+F(neighbor(ni)).x;
                            sumy=sumy+F(neighbor(ni)).y;
                            if ni~=imn
                            LL(neighbor(ni))=-mn;
                            end
                            
                        end
                        LL(neighbor(imn))=sum;
                        F(neighbor(imn)).x=sumx;
                        F(neighbor(imn)).y=sumy;
                        c(i,j)=neighbor(imn);                       
                   else
                       largest_label=largest_label+1;
                        LL(largest_label)=1;
                        c(i,j)=largest_label;
                        F(largest_label).x=j;
                        F(largest_label).y=i;
                   end
            
                    
                        
                       
                       
               
            
        end
    end
end

% h=[0; LL(1:largest_label)];
% % for i=1:largest_label
% %        
% %         c(c==i)=LL(i);
% %        
% % end
% c=h(c+1);
% labs=unique(c);
% labs=labs(2:end); %remove zero label
% csize=zeros(length(labs),1);
% for i=1:length(labs)
%     curclust=(c==labs(i));
%     csize(i)=sum(sum(curclust));
%     [y x]=find(curclust);
%     cpos(i,1)=mean(x);
%     cpos(i,2)=mean(y);
%     
% end
    LL=LL(1:largest_label);
    F=F(1:largest_label);
function r=properLabel(Sn)
    r=Sn;
    t=r;
    t=-LL(t);
    while t>0
        r=t;
        t=-LL(t);
    end
end
end
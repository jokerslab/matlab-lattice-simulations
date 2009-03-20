function [c labs csize cpos]=clusterCount(s,value)
if nargin<2
    value=-1;
end
Not=value+25;
[r,c1]=size(s);
c=zeros(size(s));
largest_label=0;
LL=zeros(r*c1,1);
for i=1:r+1
    im1=i-1;
        if im1==0,im1=1;end
        if i==r+1
            i=1;end
    for j=1:c1+1
        
        jm1=j-1;
        if jm1==0,jm1=1;end
        if j==c1+1,j=1;end
        if s(i,j)==value;
            if(im1==r)
                oldc=LL(c(i,j));
            end
            above=s(im1,j);
            if i==1 && im1==1,above=Not;end
            left=s(i,jm1);
            if j==1 && jm1==1,left=Not;end
            if (left~=value) && above~=value  && c(i,j)==0
                largest_label=largest_label+1;
                LL(largest_label)=largest_label;
                c(i,j)=largest_label;
            else
                if c(i,j)==0
                    
                if (above == value)
                    if (left==value)
                        mn=min(LL(c(im1,j)),LL(c(i,jm1)));
                        mx=max(LL(c(im1,j)),LL(c(i,jm1)));
                             %LL(mx)=mn;
                        
                        LL(LL==mx)=mn;
                        c(i,j)=LL(c(im1,j));
                       
                        %c(c==c(im1,j))=c(i,jm1);
                    else
                        c(i,j)=LL(c(im1,j));
                    end
                else
                    c(i,j)=LL(c(i,jm1));
                end
                else
                    if (above == value)
                        if (left==value)
                        mn=min([LL(c(im1,j)),LL(c(i,jm1)),LL(c(i,j))]);
                        mx=max([LL(c(im1,j)),LL(c(i,jm1)),LL(c(i,j))]);
                   
                        
                        %LL(mx)=mn;
                        
                        LL(LL==mx)=mn;

                        c(i,j)=LL(c(im1,j));
                       
                        %c(c==c(im1,j))=c(i,jm1);
                        else
                           mn=min(LL(c(im1,j)),LL(c(i,j)));
                            mx=max(LL(c(im1,j)),LL(c(i,j)));
                             %LL(mx)=mn;
                        
                        LL(LL==mx)=mn; 
                        c(i,j)=LL(c(im1,j));
                        end
                    elseif (left==value)
                        mn=min(LL(c(i,j)),LL(c(i,jm1)));
                        mx=max(LL(c(i,j)),LL(c(i,jm1)));
                             %LL(mx)=mn;
                        
                        LL(LL==mx)=mn;
                        c(i,j)=LL(c(i,jm1));
                    else
                        %Nothing
                    end
                end
            end
            if (im1==r)
                LL(LL==oldc)=LL(c(i,j));
                %c(c==oldc)=c(i,j);
            end
            
        end
    end
end

h=[0; LL(1:largest_label)];
% for i=1:largest_label
%        
%         c(c==i)=LL(i);
%        
% end
c=h(c+1);
labs=unique(c);
labs=labs(2:end); %remove zero label
csize=zeros(length(labs),1);
for i=1:length(labs)
    curclust=(c==labs(i));
    csize(i)=sum(sum(curclust));
    [y x]=find(curclust);
    cpos(i,1)=mean(x);
    cpos(i,2)=mean(y);
    
end
    

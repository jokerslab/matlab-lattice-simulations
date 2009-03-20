function [c labs csize]=clusterCount2(s)
[r,c1]=size(s);
c=zeros(size(s));
largest_label=0;
LL=zeros(r*c1,1);
for i=1:r
    im1=i-1;
        if im1==0,im1=1;end    
    for j=1:c1
        jm1=j-1;
        if jm1==0,jm1=1;end
        if s(i,j)==-1
           
            above=s(im1,j);
            if i==1 && im1==1,above=1;end
            left=s(i,jm1);
            if j==1 && jm1==1,left=1;end
            if j==c1
                right=s(i,1);
            else
                right=1;
            end
            if (left==1) && above==1 && right==1
                largest_label=largest_label+1;
                LL(largest_label)=largest_label;
                c(i,j)=largest_label;
            else
                if (above ~= 1)
                    if (left~=1)
                        
                        mn=min(LL(c(im1,j)),LL(c(i,jm1)));
                        mx=max(LL(c(im1,j)),LL(c(i,jm1)));
                        %LL(mx)=mn;
                        
                        LL(mx)=mn;
                        c(i,j)=LL(c(im1,j));
                       
                        %c(c==c(im1,j))=c(i,jm1);
                    else
                        c(i,j)=LL(c(im1,j));
                    end
                else
                    c(i,j)=LL(c(i,jm1));
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
    csize(i)=sum(sum(c==labs(i)));
end
    

function [c LL F pts]=clusterCountEHK2(s,value)
%Find the clusters of connected'values' in the given lattic state s
if nargin<2
    value=-1;
elseif length(value)>1
  for i=1:length(value)
      s(s==value(i))=value(1);
  end
  value=value(1);
end
Not=value+25;
[r,c1]=size(s);
c=zeros(size(s));
largest_label=0;
LL=zeros(r*c1,1);
periodicBCs=true;


%F=struct('x',zeros(1,r*c1),'y',zeros(1,r*c1),...
%    'x2',zeros(1,r*c1),'y2',zeros(1,r*c1));
F.x=zeros(r*c1,1);
F.y=F.x;
F.x2=F.x;
F.y2=F.x;
pts=cell(r*c1,1);
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

            
                neighbor=[0 0];
                
                if above==value
                    pl=properLabel(c(im1,j));
                   neighbor(1)=pl;
                end
                if left==value
                    neighbor(2)=properLabel(c(i,jm1));
                    
                end

               
                   if neighbor(1)>0 || neighbor(2)>0
                      
                       if neighbor(1)==neighbor(2)
                           neighbor=neighbor(1);
                       elseif neighbor(1)>0
                           if neighbor(2)>0
                               neighbor=neighbor(1:2);
                           else
                               neighbor=neighbor(1);
                           end
                       else
                           neighbor=neighbor(2);
                       end

                       %neighbor=unique(neighbor(neighbor>0));
                       
                       %[c,LL,F]=mergeClusters(i,j,neighbor,c,LL,F);
                       [mn,imn]=min(neighbor);
                        
                        sumc = 1+sum(LL(neighbor));
                        sumx = j+sum(F.x(neighbor));
                        sumy = i+sum(F.y(neighbor));
                        F.x2(neighbor(imn))=j*j+sum(F.x2(neighbor));
                        F.y2(neighbor(imn))=i*i+sum(F.y2(neighbor));
                        pts{mn}=[(j-1)*c1+i pts{neighbor}];
                        %Label of labels negative sign points to
                        %'proper label'
                        LL(neighbor(neighbor~=mn))=-mn;



                        LL(neighbor(imn))=sumc;
                        F.x(neighbor(imn))=sumx;
                        F.y(neighbor(imn))=sumy;
                        c(i,j)=neighbor(imn);            
                                  
                   else
                       largest_label=largest_label+1;
                        LL(largest_label)=1;
                        c(i,j)=largest_label;
                        F.x(largest_label)=j;
                        F.y(largest_label)=i;
                        F.x2(largest_label)=j*j;
                        F.y2(largest_label)=i*i;
                        pts{largest_label}=(j-1)*c1+i;
                   end
  
        end
    end
end
if periodicBCs
for i=1:r
    if s(i,1)==value && s(i,c1)==value
        s1=properLabel(c(i,1));
        s2=properLabel(c(i,c1));
        if s1<s2
            F.x2(s1)=F.x2(s1)+F.x2(s2)-2*c1*F.x(s2)+LL(s2)*c1*c1;
            F.y2(s1)=F.y2(s1)+F.y2(s2);
            F.y(s1)=F.y(s1)+F.y(s2);
            F.x(s1)=F.x(s1)+F.x(s2)-c1*LL(s2);
            pts{s1}=[pts{s1} pts{s2}];
            LL(s1)=LL(s1)+LL(s2);
            LL(s2)=-s1;
           
        end
        if s2<s1
            F.x2(s2)=F.x2(s1)+F.x2(s2)-2*c1*F.x(s2)+LL(s2)*c1*c1;
            F.y2(s2)=F.y2(s1)+F.y2(s2);
            F.y(s2)=F.y(s1)+F.y(s2);
            F.x(s2)=F.x(s1)+F.x(s2)-c1*LL(s2);
            pts{s2}=[pts{s1} pts{s2}];
            LL(s2)=LL(s2)+LL(s1);
            LL(s1)=-s2;
        end
        
        
    end
end

for j=1:c1
    
    if s(1,j)==value && s(r,j)==value
        s1=properLabel(c(1,j));
        s2=properLabel(c(r,j));
        if s1<s2
            F.x2(s1)=F.x2(s1)+F.x2(s2);
            F.y2(s1)=F.y2(s1)+F.y2(s2)-2*r*F.y(s2)+LL(s2)*r*r;
            F.y(s1)=F.y(s1)+F.y(s2)-r*LL(s2);
            F.x(s1)=F.x(s1)+F.x(s2);
            pts{s1}=[pts{s1} pts{s2}];
            LL(s1)=LL(s1)+LL(s2);
            LL(s2)=-s1;
        end
        if s2<s1
            F.x2(s2)=F.x2(s1)+F.x2(s2);
            F.y2(s2)=F.y2(s1)+F.y2(s2)-2*r*F.y(s2)+LL(s2)*r*r;
            F.y(s2)=F.y(s1)+F.y(s2)-r*LL(s2);
            F.x(s2)=F.x(s1)+F.x(s2);
            pts{s2}=[pts{s1} pts{s2}];
            LL(s2)=LL(s2)+LL(s1);
            LL(s1)=-s2;
        end
        
        
    end
    
end
end
    LL=LL(1:largest_label);
    pts=pts(1:largest_label);
    
    F.x=F.x(1:largest_label);
    F.y=F.y(1:largest_label);
    F.x2=F.x2(1:largest_label);
    F.y2=F.y2(1:largest_label);
    F.s=LL(1:largest_label);
    F.rs2=F.s.^-2.*((F.s.*F.x2-F.x.*F.x+F.s.*F.y2-F.y.*F.y));
    
    
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
function nc=testState(states)
n=length(states);
%s=sign(randn(40,85));

s=ones(4)*states(3);
ne=zeros(8,4,4);
    for i=1:4
        for j=1:4
            ne(:,i,j)=neighbors2d(i,j,4,4);
        end
    end
nc=[0];
for a=1:n
    for b=1:n
        for c=1:n
            for d=1:n
                %if d~=a
                    for e=1:n
                        %if e~=b
                            for f=1:n
                                %if f~=c
                                    for g=1:n
                                        for h=1:n
                                            
                                s(2,2)=states(g);
                                s(2,3)=states(h);
                                s(1,2)=states(a);
                                s(3,2)=states(b);
                                s(2,1)=states(c);                                
                                s(3,3)=states(d);
                                s(1,3)=states(e);
                                s(2,4)=states(f);
                                p=unique(classifypairs2(s,ne));
                                        for i=1:length(p)
                                        if sum(find(nc==p(i)))==0
                                            nc=[nc; p(i)];
                                        end 
                                        end 
                                        
                                        
                                        end
                                        
                                    end
                                %end
                            end
                        %end
                    end
                %end
            end
        end
    end
end
nc=unique(nc);
%   
s=zeros(3,4*3^8);
s(2:12*3:12*3^8)=1;
s(2+12:12*3:12*3^8)=2;
s(2+24:12*3:12*3^8)=3;
s(4:12:3
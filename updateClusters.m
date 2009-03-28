%Update Labeled Cluster Matrix
if cs(si)>0  %site occupied before swap
   preCSite=si;
   postCsite=n_1(nd);

else
  preCSite=n_1(nd);
  postCsite=si;
end
%count the neighbors before swap
   nbs=sum(cs(neighbors(1:4,preCSite))>0);
    if nbs>1
            dan=neighbors(1:4,preCSite);
            dan=dan(cs(dan)>0);
    end
   %Swap
   temp=cs(si);
   cs(si)=cs(n_1(nd));
   cs(n_1(nd))=temp;

   %get neighbors after swap
   nas=properLabel(cs(neighbors(1:4,postCsite)),LL);
   cspl=properLabel(cs(postCsite),LL);
    [i,j]=ind2sub(size(cs),postCsite);
   if nbs==0
       if any(nas>0)
           nas=unique(nas(nas>0));

           %merge clusters
           [mn,imn]=min(nas);
        
            LL(nas(imn)) = 1+sum(LL(nas));
            F.x(nas(imn)) = j+sum(F.x(nas));
            F.y(nas(imn)) = i+sum(F.y(nas));
            F.x2(nas(imn))=j*j+sum(F.x2(nas));
            F.y2(nas(imn))=i*i+sum(F.y2(nas));

            %Label of labels negative sign points to
            %'proper label'
            LL(nas(nas~=mn))=-mn;
            LL(cspl)=0;
            cs(postCsite)=nas(imn);            
            
            %can reuse the old label
            %need to write a function to handle this
            
       else
           % nothing happens
       end
   else
       if any(nas==cspl)
           if any(nas(nas~=cspl)>0)
               %merge new cluster into us
               nas=unique(nas(nas~=cspl & nas>0));
               LL(cspl)=LL(cspl)+sum(LL(nas));
               
               LL(nas)=-cspl;
           end
       else
           
           %substract moving piece from old cluster
           LL(cspl)=LL(cspl)-1;
           
           if any(nas>0)
               nas=unique(nas(nas>0));

           %merge clusters
               [mn,imn]=min(nas);
               
                LL(nas(imn)) = 1+sum(LL(nas));
                F.x(nas(imn)) = j+sum(F.x(nas));
                F.y(nas(imn)) = i+sum(F.y(nas));
                F.x2(nas(imn))=j*j+sum(F.x2(nas));
                F.y2(nas(imn))=i*i+sum(F.y2(nas));

                %Label of labels negative sign points to
                %'proper label'
                LL(nas(nas~=mn))=-mn;

                cs(postCsite)=nas(imn);     
           else
               %we broke off our old cluster and formed newone
               largest_label=largest_label+1;
                LL(largest_label)=1;
                cs(postCsite)=largest_label;
                F.x(largest_label)=j;
                F.y(largest_label)=i;
                F.x2(largest_label)=j*j;
                F.y2(largest_label)=i*i;
%                pts{largest_label}=postCsite;
               
               
           end
           
       end
        if nbs>1
            
            line=scanBuild([],dan(1),neighbors(1:4,:),s);
           
            if length(dan)>2
                %debugme
                dan2=dan(~ismember(dan,line));
                if ~isempty(dan2)
                    largest_label=largest_label+1;
                    cs(line)=largest_label;
                        LL(largest_label)=length(line);
                        LL(cspl)=LL(cspl)-length(line);
                        [i,j]=ind2sub(size(cs),line);
                        F.x(largest_label)=sum(j);
                        F.y(largest_label)=sum(i);
                        F.x2(largest_label)=sum(j.*j);
                        F.y2(largest_label)=sum(i.*i);
                        
                        if length(dan2)>1
                            line=scanBuild([],dan2(1),neighbors(1:4,:),s);
                            dan=dan2;
                        else
                            line=dan(2);
                        end
                end
             end
                 if any(dan(2)==line)
                %ok both points still in cluster
                 else
                        largest_label=largest_label+1;
                        cs(line)=largest_label;
                        LL(largest_label)=length(line);
                        LL(cspl)=LL(cspl)-length(line);
                        [i,j]=ind2sub(size(cs),line);
                        F.x(largest_label)=sum(j);
                        F.y(largest_label)=sum(i);
                        F.x2(largest_label)=sum(j.*j);
                        F.y2(largest_label)=sum(i.*i);
                 end
            
            
            
        end
   end

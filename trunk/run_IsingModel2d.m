count=1;
for i=1:5
    kbT=10^(i/2-1.5);
    k=0;
    for phi=.95:-.05:.5
        k=k+1;
        for repeat=1:10
             [s{count} Energy{count} cor{count} clust{count} cave{count}]= ...
             Ising2D(kbT,phi);
            [a b c]=clusterCount(s{count});
            m1(repeat)=max(c);
            nc1(repeat)=length(c);
            mn1(repeat)=mean(c);
%            if repeat==3
%                [a b c]=clusterCount(s{count});
%                 plotLattice(a,b)
%                 title(['Position after 20000 steps. T=' num2str(kbT) ' \phi=' num2str(1-phi)])
%                 pause
%            end
               count=count+1;
             pause(1)
             close all

%Plot the energies after loading the data
            %[e1 e2 e3 e4 e5 e6 e7 e8 e9 e10]=Energy{1:10};
            %t=0:20000;
            %plot(t,e1,t,e2,t,e3,t,e4,t,e5,t,e6,t,e7,t,e8,t,e9,t,e10)
        end
       mnBigC(i,k)=mean(m1);
       sdBigC(i,k)=std(m1);
       mn(i,k)=mean(mn1);
       sdmn(i,k)=std(mn1);
       nc(i,k)=mean(nc1);
       sdnc(i,k)=std(nc1);
        
    end
    i
end
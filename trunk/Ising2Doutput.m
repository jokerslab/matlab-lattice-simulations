if saveFrames
    filename = [directory, 'data.mat'];
    save(filename);
end

figure()
plotLattice(s, [1 2 3],{'none' [.5 .5 .5],'black'})
title(['State after ' num2str(t) ' steps'])
if saveplot
    
    print('-r600','-depsc',[filestr 'Finalpos']);
end
%figure()
%plotLattice(pairs,[0 1 2 3 4 5 6 -1])
figure()
%plotLattice(classifypairs(s),[0 1 2 3 4 5 6 -1])
plot(0:t,Energy(1:t+1));
ylabel('Energy')
xlabel('Steps')
title('Energy of System')
if saveplot
filestr=['nucMC/pics/' startTime 'Energy'];
print('-r600','-depsc',filestr);
end
ts=0:Sps:t;
if checkAbsorb
figure
plot(relt,release)
end
%figure
%plot(ts,cor(5,6,:))
%title('Correlation 1 lattice site away')
%filestr=['nucMC/pics/' startTime 'Correlation'];
%print('-r600','-depsc',filestr);
%figure
%plot(ts,cave)
%title('Average cluster (>10) size')
%if saveplot
%filestr=['nucMC/pics/' startTime 'Csize'];
%print('-r600','-depsc',filestr);
%end

% function pausefig(src,evnt)
%       if evnt.Character == 't'
%          pause(1)
%       end
% end
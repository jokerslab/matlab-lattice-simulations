%Produce the Graphs and save the data from a LatticeSimulation
if params.saveFrames
    filename = [directory, 'data.mat'];
    save(filename);
end

figure()
plotLattice(s, [1 2 3],{'none' [.5 .5 .5],'black'})
title(['State after ' num2str(t) ' steps'])
if params.saveplot
    
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
if params.saveplot
filestr=['nucMC/pics/' startTime 'Energy'];
print('-r600','-depsc',filestr);
end
ts=0:params.Sps:t;
if checkAbsorb
figure
plot(relt,release)
end
figure();
plot(ts,fre1/4096,'b',1:t,frac3/4096,'r',1:t,frac2/4096,'g')
xlabel('Simulation Steps');
ylabel('Concentration')
legend('Free B','B*','B')
%title('k_bT = 1 \chi_{AB} = -0.1 \chi_{AB^*} = 0.2 \epsilon = -0.5 \gamma= 0.5')
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
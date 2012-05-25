function [ et ] = computeInteractionParameters( chiSA,chiSC,epsilon,gamma )
%Calcuate the individual pair interactions between states given simulation
%parameters
% Output pairwise energy components in a matrix (ie addressable by state)
%Given each state as an integer (1=Solvent, 2=Amorphous, 3=Crystal,
%4=inert, 5=Water...) the energy transfer matrix (et) is a symmetric matrix
%where the energy between two states is given by the value indexed in the
%matrix from state i, and state j. The et(state_i,state_b) returns the w_ij
%interaction energy between state i and j
% The et matrix contains more parameters than necessary for the system, the
% system is defined by several interaction parameters and these are used to
% define the relative energies between states to make the simulation
% computation easier.
% Given the interaction parameters,surface energy and crystal bulk
%   compute the relative indivdual interaction energies
SS=0; 
AA=0;
SA=chiSA;
SC=gamma;

AC=chiSC+epsilon/2;
CC=epsilon;

WW=0;SW=0;CW=0;AW=0;  %Wall interactions

%Water Interactions

AQS=0;AQA=.6;AQC=.1;AQW=0;AQAQ=0;

et=-1*[SS SA SC SW AQS;
    SA AA AC AW AQA;
    SC AC CC CW AQC;
    SW AW CW WW AQW;
    AQS AQA AQC AQW AQAQ];
end


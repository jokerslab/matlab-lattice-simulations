%Bond energies
%Parameter file that sets the Energy between different neighbor states
%Energy interactions stored as a symmetric matrix for lookup
SS=0;  %1
SA=.1;  %2 chiSA=-.1
SC=-1.1; %-2 gamma=1.1
AA=0; %1
AC=-.1; %-1 chicAC=.6
CC=1; %4 epsislon=-1

SS=-1;  %1
SA=-.9;  %2
SC=-2.1; %-2
AA=-1; %1
AC=-1.1; %-1
CC=0; %4


WW=0;SW=0;CW=0;AW=0;  %Wall interactions

%Water Interactions

AQS=0;AQA=.6;AQC=.1;AQW=0;AQAQ=0;

et=[SS SA SC SW AQS;
    SA AA AC AW AQA;
    SC AC CC CW AQC;
    SW AW CW WW AQW;
    AQS AQA AQC AQW AQAQ];
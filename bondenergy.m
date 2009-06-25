%Bond energies
SS=0;  %1
SA=.1;  %2
SC=-1.1; %-2
AA=0; %1
AC=-.1; %-1
CC=1; %4

WW=0;SW=0;CW=0;AW=0;  %Wall interactions

%Water Interactions

AQS=0;AQA=.6;AQC=.1;AQW=0;AQAQ=0;

et=[SS SA SC SW AQS;
    SA AA AC AW AQA;
    SC AC CC CW AQC;
    SW AW CW WW AQW;
    AQS AQA AQC AQW AQAQ];
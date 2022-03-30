%% Junction Potential
cSAMPLE = {'Hilus', 'Cortex'};

for iSmp = 1:length(cSAMPLE)
    TC = 24;
    R=8.314;
    F=96485;
    TK=273.15 + TC;
    
    sIon=struct;
    
    sIon(1).name='Sodium';
    sIon(1).z= 1;%Valence
    sIon(1).u= 0.682; %mobility
    sIon(1).CPipette= 0;%intracellular concentration (unit must match COut)
    %extracellular concentration (unit must match CIn)
    if iSmp == 1, sIon(1).CACSF=151.2; else, sIon(1).CACSF=153.25; end
    
    sIon(2).name='Potassium';
    sIon(2).z= 1;%Valence
    sIon(2).u= 1; %mobility
    sIon(2).CPipette= 144;%intracellular concentration (unit must match COut)
    sIon(2).CACSF= 2.5;%extracellular concentration (unit must match CIn)
    
    sIon(3).name='Chloride';
    sIon(3).z= -1;%Valence
    sIon(3).u= 1.0388; %mobility
    sIon(3).CPipette= 6;%intracellular concentration (unit must match COut)
    %extracellular concentration (unit must match CIn)
    if iSmp == 1, sIon(3).CACSF= 134.5; else, sIon(3).CACSF= 136.5; end
    
    sIon(4).name='Magnesium';
    sIon(4).z= 2;%Valence
    sIon(4).u= 0.361; %mobility
    sIon(4).CPipette= 3;%intracellular concentration (unit must match COut)
    %extracellular concentration (unit must match CIn)
    if iSmp == 1, sIon(4).CACSF= 3; else, sIon(4).CACSF= 1; end
    
    sIon(5).name='Calcium';
    sIon(5).z= 2;%Valence
    sIon(5).u= 0.4048; %mobility
    sIon(5).CPipette= 0;%intracellular concentration (unit must match COut)
    sIon(5).CACSF= 2;%extracellular concentration (unit must match CIn)
    
    sIon(6).name='Gluconate';
    sIon(6).z= -1;%Valence
    sIon(6).u= 0.33; %mobility
    sIon(6).CPipette= 144;%intracellular concentration (unit must match COut)
    sIon(6).CACSF= 0;%extracellular concentration (unit must match CIn)
    
    sIon(7).name='HEPES';
    sIon(7).z= -1;%Valence
    sIon(7).u= 0.3; %mobility
    sIon(7).CPipette= 10;%intracellular concentration (unit must match COut)
    sIon(7).CACSF= 0;%extracellular concentration (unit must match CIn)
    
    sIon(8).name='HCO3';
    sIon(8).z= -1;%Valence
    sIon(8).u= 0.605; %mobility
    sIon(8).CPipette= 0;%intracellular concentration (unit must match COut)
    %extracellular concentration (unit must match CIn)
    if iSmp == 1, sIon(8).CACSF= 26; else, sIon(8).CACSF= 26.2; end
    
    sIon(9).name='H2PO4';
    sIon(9).z= -1;%Valence
    sIon(9).u= 0.45; %mobility
    sIon(9).CPipette= 0;%intracellular concentration (unit must match COut)
    %extracellular concentration (unit must match CIn)
    if iSmp == 1, sIon(9).CACSF= 1; else, sIon(9).CACSF= 1.25; end
    
    SFup=0; SFdown=0; SP=0; SS=0;
    for i=1:length(sIon)
        SFup=SFup+(sIon(i).z*sIon(i).u*(sIon(i).CACSF-sIon(i).CPipette));
        SFdown=SFdown+((sIon(i).z^2)*sIon(i).u*(sIon(i).CACSF-sIon(i).CPipette));
        SP=SP+((sIon(i).z^2)*sIon(i).u*(sIon(i).CPipette));
        SS=SS+((sIon(i).z^2)*sIon(i).u*(sIon(i).CACSF));
    end
    SF=SFup/SFdown;
    
    VJunc=1000*((R*TK)/F)*SF*log(SP/SS);
    fprintf('%s:\tJunction potential: %s mV\n', cSAMPLE{iSmp}, num2str(VJunc,3))
end
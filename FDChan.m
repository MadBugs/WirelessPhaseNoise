classdef FDChan < matlab.System
    % Frequency-domain multipath channel
    properties
        % Configuration
        carrierConfig;   % Carrier configuration
        waveformConfig;  % Waveform configuration
        
        % Path parameters
        gain;  % Path gain in dB
        dly;   % Delay of each path in seconds
        aoaAz, aoaEl; % Angle of arrival of each path in degrees      
         fd;    % Doppler shift for each path
        
        rxVel = [30,0,0]';  % Mobile velocity vector in m/s
        fc = 28e9;    % Carrier freq in Hz
        
         gainComplex;  % Complex gain of each path
        
        % SNR parameters
        Etx = 1;       % average energy per PDSCH symbol 
        EsN0Avg = 20;  % Avg SNR per RX symbol in dB
       
        % Symbol times
         symStart;  % symStart(i) = start of symbol i relative to subframe      
         
        % Verbose
        Verb=1; 
    end
    methods
        function obj = FDChan(carrierConfig, varargin)
            % Constructor
            
            % Save the carrier configuration
            obj.carrierConfig = carrierConfig;
                               
            % Set parameters from constructor arguments
            if nargin >= 1
                obj.set(varargin{:});
            end
            
            % ***TODO:  Create complex path gain for each path
            %   obj.gainComplex = ... 
            % obj.gain, obj.dly, obj.fc
            len = length(obj.gain);
            obj.gainComplex = db2mag(obj.gain) .* exp(1i* -2*pi*rand(len,1));
                  
            % Compute unit vector in direction of each path
            az=obj.aoaAz'; el=obj.aoaEl';
            [ux uy uz]=sph2cart(deg2rad(az),deg2rad(el),ones(1,length(az)));
            U = [ux; uy; uz]; %U, diag(U'*U) %Check
            
            % ***TODO:  Compute the Doppler shift for each path
            %    obj.fd = ...          
            vc     = physconst('lightspeed');
            obj.fd = ((obj.rxVel'*U)*(-obj.fc/vc))';  
            %obj.fd = ((obj.rxVel'*U)*( obj.fc/vc))'; 
            
            % ***TODO:  Compute the vector of 
            % symbol times relative to the start of the subframe
            %    obj.symStart = ...  
            
            % 1) These are the relatives times to END of each symbol
            aux =cumsum(obj.waveformConfig.SymbolLengths/ ...
                        obj.waveformConfig.SampleRate);
            
            % 2) Now we R-SHIT to get the time to BEGINING of the subframe
            obj.symStart = [0 aux(1:end-1)]; 
            
            if (obj.Verb) fprintf("FDChan->FDChan() Called"); end
        end 
    end
    
    methods (Access = protected)
        
        function [rxGrid, chanGrid, noiseVar] = ...
                               stepImpl(obj, txGrid, sfNum, slotNum)
            
            % Applies a frequency domain channel and noise
            %
            % Given the TX grid of OFDM REs, txGrid, the function
            % *  Computes the channel grid, chanGrid, given the 
            %    subframe number <sfNum> and slotNum <slotNum>
            %
            % *  Computes the noise variance per symbol, noiseVar,
            %    for a target SNR
            %
            % *  Applies the channel and noise to create the RX grid 
            %    of symbols, rxGrid.
    
            % ???TODO           
            E2_chan  = obj.gainComplex' * obj.gainComplex; % <G*,G>
            noiseVar = obj.Etx * E2_chan / db2pow(obj.EsN0Avg);
                     
            JJ = length(txGrid(1,:)); 
            NN = length(txGrid(:,1));
            jj = [0:JJ-1]';  
            nn = [0:NN-1]';     
                     
            % Each <sfNum>   range:[0:10-1] adds (sfNum * 1ms)  of delay  
            % Each <slotNum> range:[0:SlotsPerSubframe-1] adds 
            %   (symStart(slotNum*SlotsPerSubframe+1) * 1ms)    of delay
            tsymbs = sfNum * 1e-3     + ...
                     obj.symStart(jj  + ...
                         slotNum * obj.waveformConfig.SymbolsPerSlot+1);
            
            scs      = obj.carrierConfig.SubcarrierSpacing*1e3;           
            dly_term = (nn*scs .* obj.dly');
  
            % Now we must create the chanGrid
            chanGrid=zeros(NN, JJ);
            for col_j  = [1:JJ]
              gains_kj = obj.gainComplex' .* ...  % will use bcast below
                         exp(2*pi*1i*(tsymbs(col_j)*obj.fd' + dly_term));
              chan_j   = sum(gains_kj, 2); % sums over k cols
              chanGrid(:, col_j) = chan_j;
            end
            
            % Now creates the noise and creates the rxGrid
            Bruit = (randn(NN,JJ) + 1i*randn(NN,JJ)) * sqrt(noiseVar/2);          
            rxGrid= chanGrid.*txGrid + Bruit;  
            
            %fprintf("FDChan->stepImpl() Called");

%           UNIT TEST: Looks right. Tested for
%           sfNum  =9; %Range: [0:10-1]
%           slotNum=7; %Range: [0: 8-1]
%           And got almost a SUBFRAME TIME of 1ms AS EXPECTED!
%
%             format long; 
%             tsymbs'
% 
%             ans =            
%                0.009875130208333
%                0.009884049479167
%                0.009892968750000
%                0.009901888020833
%                0.009910807291667
%                0.009919726562500
%                0.009928645833333
%                0.009937565104167
%                0.009946484375000
%                0.009955403645833
%                0.009964322916667
%                0.009973242187500
%                0.009982161458333
%                0.009991080729167

        end   
     end
end
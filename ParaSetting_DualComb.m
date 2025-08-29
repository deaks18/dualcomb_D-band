FLAG.if_PlotFigures = 0;
FLAG.if_PlotMIMO = 0;

FLAG.if_CreateDACFiles = 0;
FLAG.if_CaptureNew = 1;
FLAG.if_PlotSYNC = 1;
FLAG.Tx_Conjugation = 0;
FLAG.if_RunSymbolDecision = 1;
FLAG.if_saveMIMOfilter = 0;
FLAG.if_PlotConstellation = 1;

FLAG.N_CAPTURE_START = 1;

EQUIPMENT.SaDAC = 2.4e9; %120e9; %88e9;
EQUIPMENT.SoftSa = EQUIPMENT.SaDAC;
EQUIPMENT.Fs = 10e9;

SC.RRC_Shaping = 0.1;
SC.SymbolRates = 200e6;
SC.Mode = 'QAM';
SC.ModulationFormats = 64; % M-QAM
SC.DSP.NumTaps = 201;

SC.M = SC.ModulationFormats;
SC.Ini_Phase = 0;
SC.band_center_freq = 0; %0; %0;
SC.Symbol_Order = 'gray';
SC.DimensionMode = 'PS-QAM';
SC.Betas = 1;
SC.beta_actual = 1;

SC.DSP.mu = 6e-4;
SC.DSP.TrainingSymbols = 1000;
SC.DSP.PLLg = 0.001;
SC.DSP.N_equaliser_loops = 4;

SC.N_padding = 0;
SC.real_MIMO = 0;

FLAG.if_CalculateMultiframes = 0;
FLAG.N_POL = 1;
EQUIPMENT.N_BANDS = 1;
SC.N_bands = 1;

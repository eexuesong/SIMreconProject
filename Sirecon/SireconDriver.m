%{
All calculation starts with a SIM_Reconstructor object

As seen in cudaSireconDriver.cpp, SIM reconstruction is done in the following steps:
1. Instantiate a SIM_Reconstructor object by passing along the command line;
2. For all the time points, do:
2a. Load raw data and preprocess
2b. call member processOneVolume() to reconstruct the current time point
2c. write reconstruction result of current time point
3. Close file(s) and exit
%}

% g = gpuDevice(1); reset(g);
% disp(['GPU Memory: ',num2str(g.FreeMemory / 1024 / 1024 / 1024), ' GB']);

clc;
clear all;

arg = struct(...
    'ifiles', 'D:\SIMreconProject\Sirecon\DataForTest',...      % input file (or data folder in TIFF mode)
    'ofiles', 'U2OS_Tomm20_3DSIM(1P35NA)',...                       % output file (or filename pattern in TIFF mode)
    'otffiles', 'D:\SIMreconProject\Sirecon\OTF\OTF_3DSIM_1P35NA_488.tif',... % OTF file    
    'usecorr', [],...                   % flat-field correction file provided
    'ndirs', 3,...                      % number of directions
    'nphases', 5,...                    % number of phases per direction
    'nordersout', 0,...                 % number of output orders; must be <= norders
    'angle0', -0.19142,...              % angle of the first direction in radians (1.27 NA: -1.3994; 1.35 NA:-0.19142) 
    'ls', [0.1993, 0.1994, 0.1994],...  % line spacing of SIM pattern in microns
    'na', 1.35,...                      % Detection numerical aperture
    'nimm', 1.406,...                   % refractive index of immersion medium
    'nmed', 1.406,...                   % refractive index of sample medium
    'iter_input', 1,...                 % iteration number for input RL deconvolution
    'iter_output', 2,...                % iteration number for output RL deconvolution
    'zoomfact', 2,...                   % lateral zoom factor
    'explodefact', 1.0,...              % artificially exploding the reciprocal-space distance between orders by this factor
    'zzoom', 1,...                      % axial zoom factor
    'background', 100.0,...             % camera readout background (be will be subtracted from input images)
    'wiener', 0.001,...                 % Wiener constant
    'linear_wiener', [],...             % Linear increment of Wiener constant for time-lapse imaging
    'forcemodamp', [],...               % modamps forced to these values
    'k0angles', [2.948954, 1.900853, -2.288121],... % user given pattern vector k0 angles for all directions
    'gammaApo', 1.0,...                   % output apodization gamma; 1.0 corresponds to a triangular shape, lower numbers (e.g. 0.5) decrease medium frequency response (roughly as A(k) = 1 - k ^ b with k = 0...1)
    'saveprefiltered', [],...           % save separated bands (in frequency domain) into a file and exit
    'savealignedraw', [],...            % save drift-fixed raw data (half Fourier space) into a file and exit
    'saveoverlaps', [],...              % save overlap0 and overlap1 (real-space complex data) into a file and exit
    'config', [],...                    % name of a file of a configuration.
    'lambda', 488,...                   % excitation wavelength: 488, 561
    'wavelength', 525);                 % emission wavelength: 525, 607 (only used for TIFF files)


% 1.27 NA, 488 nm: 'ls', [0.20909913, 0.20898064, 0.2090098]
% 1.27 NA, 561 nm: 'ls', [0.24000179, 0.24009595, 0.24017175]
% 1.27 NA:         'k0angles', [1.74091516, 2.78794151, 0.68857549]

% 1.35 NA, 488 nm: [0.1993, 0.1994,0.1994]
% 1.35 NA, 561 nm: [0.2287, 0.2287, 0.2287]
% 1.27 NA:         'k0angles', [2.948954, 1.900853, -2.288121]


% arg.('useRLonInput') = [];               % use 2D Richardson-Lucy deconvolution on input images
% arg.('useRLonOutput') = [];              % use 2D / 3D Richardson-Lucy deconvolution on output images
% arg.('nofilteroverlaps') = [];          % do not filter the overlaping region between bands, usually used in trouble shooting
arg.('otfRA') = [];                     % using rotationally averaged OTF
arg.('otfPerAngle') = [];               % using one OTF per SIM angle
arg.('fastSI') = [];                    % SIM data is organized in Z->Angle->Phase order; default being Angle->Z->Phase
arg.('searchforvector') = [];           % search for k0 at the 0th time point (If warning message appears, comment this property)
% arg.('k0searchAll') = [];               % search for k0 at all time points
arg.('equalizez') = [];                 % bleach correcting for z
arg.('equalizet') = [];                 % bleach correcting for time
arg.('dampenOrder0') = [];              % dampen order-0 in final assembly
% arg.('nosuppress') = [];                % do not suppress DC singularity in final assembly (good idea for 2D/TIRF data)
arg.('dampenFactor') = 1000;            % When dampening DC singularity, suppress how many times?
% arg.('nokz0') = [];                     % do not use kz=0 plane of the 0th order in modamp fit and assemblerealspace()
% arg.('twolenses') = [];                 % 4-beam SIM data
% arg.('applyOtfBeforeShift') = [];       % Multiply bands with OTF before frequency shift
% arg.('keepnegativevalues') = [];        % When saving wiener shifted images, keep negative values or not
% arg.('savefft') = [];                   % Whether to save fft data 
arg.('fftkeepnegativevalues') = [];     % When saving fft data, based on negative values or not
% arg.('GPUcomputing') = [];            % GPU acceleration or not

diary(strcat(arg.('ifiles'), "\info_488_3D.txt"));
diary on;

% Instantiate a SIM_Reconstructor object by passing along the command line
myreconstructor = SIM_Reconstructor(arg);

ls = zeros(myreconstructor.getNTimes(), 3);
k0angles = zeros(myreconstructor.getNTimes(), 3);

for it = 0 : (myreconstructor.getNTimes() - 1)
    for iw = 0 : 0
        tic;
        myreconstructor = myreconstructor.loadAndRescaleImage(it, iw);
        
        myreconstructor = myreconstructor.setCurTimeIdx(it);
        
        %{
        The main processing occurs inside this function:
		  1. Re-initialize all key device buffers under m_reconData;
          2. Fine-tune k0 vectors, modulation amplitudes and phases for all directions and orders
          3. For all directions, pre-filter separated bands, inverse FFT, and assemble the bands
        %}
        myreconstructor = myreconstructor.processOneVolume(it);
        
        myreconstructor = myreconstructor.writeResult(it, iw);
        
        myreconstructor.m_myParams.wiener = myreconstructor.m_myParams.wiener + myreconstructor.m_myParams.linearWiener;
        toc;
        
    end
    
    for direction = 1:arg.('ndirs')
        k0 = myreconstructor.m_reconData.k0(direction);
        k0mag = sqrt(k0.x * k0.x + k0.y * k0.y);
        k0angle = atan2(k0.y, k0.x);
        ls(it + 1, direction) = (myreconstructor.m_imgParams.ny * myreconstructor.m_imgParams.dy) / k0mag / 2;
        k0angles(it + 1, direction) = k0angle;
    end
end

diary off;


   


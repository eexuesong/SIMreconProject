classdef SIM_Reconstructor
    properties
        m_all_matching_files = {};
        
        m_config_file = [];
        
        m_myParams = struct(...
            'k0startangle', [],...          % Angle of the first direction in radians, default = -0.008 rad
            'linespacing', [],...           % Line spacing of SIM pattern in microns, default = 0.198802 um
            'na', [],...                    % Detection numerical aperture, default = 1.35
            'nimm', [],...                  % Refractive index of immersion medium, default = 1.406
            'nmed', [],...                  % Refractive index of sample medium, default = 1.333
            'ndirs', [],...                 % Number of directions, default = 3
            'nphases', [],...               % Number of phases per direction, default = 5
            'norders_output', [],...        % Number of output orders; must be <= norders
            'norders', [],...
            'phaseSteps', [],...            % User-specified non-ideal phase steps (i.e., not 2pi/5), one for each orientation, default = 0
            'bTwolens', [],...              % Whether to process I{^5}S dataset, default = 0
            'bFastSIM', [],...              % Fast SIM data is organized differently. SIM data is organized in Z->Angle->Phase order; default being Angle->Z->Phase
            'bRLonInput', [],...            % Wheter to use 2D Richardson-Lucy deconvolution on input images
            'iter_input', [],...            % Iteration number for input RL deconvolution
            'bRLonOutput', [],...           % Wheter to use 2D / 3D Richardson-Lucy deconvolution on output images
            'iter_output', [],...           % Iteration number for output RL deconvolution
            'zoomfact', [],...              % Lateral zoom factor, default = 2
            'z_zoom', [],...                % Axial zoom factor, default = 1
            'nzPadTo', [],...               % Pad zero sections to this number of z sections, default = 0
            'bPadBlur', [],...              % Whether to pad using 2D blurred images, default = false
            'explodefact', [],...           % Artificially exploding the reciprocal-space distance between orders by this factor, default = 1.0
            'bFilteroverlaps', [],...       % default = 1; do not filter the overlaping region between bands usually used in trouble shooting
            'recalcarrays', [],...          % Whether to calculate the overlaping regions between bands just once or always; used in fitk0andmodamps(), default = 1
            'napodize', [],...              % Whether and how many pixels to apodize (or edge softening) every 2D slice. Default = 10
            'forceamp', [],...
            'k0angles', [],...
            'bSearchforvector', [],...      % Whether to find initial estimate of modulation wave vector m_reconData.k0(direction) by cross-correlation, or just assume k0 vector known, so just fit for the modulation amplitude and phase; default = 1
            'bUseTime0k0', [],...           % Whether to use time 0's k0 fit for the rest of a time series. Search for k0 at all time points? Default = 1
            'apodizeoutput', [],...         % 0-no apodize; 1-cosine apodize; 2-triangle apodize; used in filterbands()
            'apoGamma', [],...              % Output apodization gamma; 1.0 means triangular apo, default = 1
            'bSuppress_singularities', [],...% Whether to dampen the OTF values near band centers, default = 1; used in filterbands(). Do not suppress DC singularity in final assembly (good idea for 2D/TIRF data)
            'dampenFactor', [],...          % When dampening DC singularity, suppress how many times
            'suppression_radius', [],...    % If suppress_singularities is 1, the range within which suppresion is applied, default = 10; used in filterbands()
            'bDampenOrder0', [],...         % Whether to dampen the OTF values near band centers, default = 0; used in filterbands(). Dampen order-0 in final assembly
            'bFitallphases', [],...         % In NL SIM, whether to use fitted phase for all orders or to infer higher order phase from order 1 phase fil, default = 1
            'do_rescale', [],...            % Fading correction method: 0-no correction; 1-with correction. Default = 1
            'equalizez', [],...             % Bleach correcting for z, default = 0
            'equalizet', [],...             % Bleach correcting for time, default = 0
            'bNoKz0', [],...                % default = 0; if true, no kz=0 plane is used in modamp fit and assemblerealspace()
            'botfBeforeShift', [], ...      % Whether to multiply bands with OTF before frequency shift
            'wiener', [],...                % Wiener constant, default = 0.001
            'linearWiener', [],...          % Linear increment of Wiener constant for time-lapse imaging
            'bUseEstimatedWiener', [],...   % if (m_myParams.wiener > 0) m_myParams.bUseEstimatedWiener = false
            'bRadAvgOTF', [],...            % Is radially-averaged OTF used? Default = 0
            'bOneOTFperAngle', [],...       % One OTF per SIM angle (instead of common OTF for all angles)?, default = 0
            'bFixdrift', [],...             % Whether nor not to correct drift between pattern directions, default = 0
            'drift_filter_fact', [],...     % Fraction of the NA; used as the cutoff frequency in high-pass filtering in drift estimation
            'constbkgd', [],...             % Camera readout background, default = 0.0
            'bBgInExtHdr', [],...           % In Andor EMCCD, background varies with each exposure, esp. in EM mode. Hidden-behind-aluminum-foil pixels can be used to estimate background of each exposure and stored in the extended header. When this option is true, the 3rd float of the extended header stores such estimated background values
            'bUsecorr', [],...              % Whether to use a camera flat-fielding (or correction) file, default = 0
            'corrfiles', [],...             % Name of the camera correction file if bUsecorr is 1
            'bMakemodel', [],...            % Whether to fake an ideal point source and obtain a recontruction of it (i.e., an effective PSF in some sense), default = 0
            'bSaveSeparated', [],...        % Whether to save separated bands and then quit before filterbands, default = 0
            'fileSeparated', [],...
            'bSaveAlignedRaw', [],...       % Whether to save drift-corrected raw images (within each direction) and then quit before filterbands
            'fileRawAligned', [],...
            'bSaveOverlaps', [],...         % Whether to save makeoverlaps() output into a file
            'fileOverlaps', [],...
            'ifiles', [],...
            'ofiles', [],...
            'otffiles', [],...
            'lambda', [],...
            'bkeepnegative', [],...         % Whether to keep negative values when saving wiener shifted images
            'bsavefft', [],...              % Whether to save fft data
            'bfftkeepnegative', [])         % When saving fft data, based on negative values or not
        
        m_imgParams = struct(...
            'nx_input', [],...              % m_imgParams.nx_input = header.nx; x pixel number of input image
            'ny_input', [],...              % m_imgParams.ny_input = header.ny; y pixel number of input image
            'nx', [],...                    % x pixel number of resized image padded with 0
            'ny', [],...                    % y pixel number of resized image padded with 0
            'nz', [],...                    % m_imgParams.nz /= params.nphases * params.ndirs; z pixel number of input image
            'nz0', [],...                   % pad zero sections to this number of z sections
            'nwaves', [],...                % m_imgParams.nwaves = header.num_waves; number of how many wavelengths were acquired
            'wave', [],...                  % emission wavelength
            'ntimes', [],...                % How many time points have been recorded for m_myParams.ifiles
            'curTimeIdx', [],...
            'dx', [],...                    % x direction pixel size (in microns)
            'dy', [],...                    % y direction pixel size (in microns)
            'dz', [],...                    % z direction pixel size (in microns)
            'dkr', [],...                   % Radial-direction inverse microns per pixel in data
            'dkz', [],...                   % Axial-direction inverse microns per pixel in data
            'krscale', [],...               % Ratio of radial direction pixel scales of data and otf
            'kzscale', [],...               % Ratio of axial direction pixel scales of data and otf
            'rdistcutoff', [],...           % OTF support radial limit in data pixels 
            'zdistcutoff', [],...           % OTF support axial limit in data pixels
            'apocutoff', [],...             % OTF support radial limit after reconstruction
            'zapocutoff', [],...            % OTF support axial limit after reconstruction
            'inscale', []...                % inscale is used to normalize all image layers so that their integral is 1
            );
        
        m_reconData = struct(...
            'otf', [],...
            'background', [],...
            'slope', [],...
            'backgroundExtra', [],...
            'savedBands', [],...
            'sepMatrix', [],...
            'noiseVarFactors', [],...
            'overlap0', [],...
            'overlap1', [],...
            'k0', [],...
            'k0_time0', [],...
            'k0guess', [],...
            'amp', [],...
            'sum_dir0_phase0', [],...
            'bigbuffer', [],...
            'outbuffer', [],...
            'SIM_OTF', [],...
            'RLbuffer', []...
            );
        
        m_driftParams = struct(...
            'driftlist', struct('x', 0, 'y', 0, 'z', 0),...
            'phaseList', 0,...
            'driftAcq', struct('x', 0, 'y', 0),...
            'drifts', struct('x', 0, 'y', 0, 'z', 0),...
            'phaseAbs', 0,...
            'timestamps', 0,...
            'timestamps_for_fitting', 0,...
            'expDose', 0,...
            'drift_bt_dirs', struct('x', 0, 'y', 0, 'z', 0));
        
        otf = [];        
        m_zoffset = [];
    end
    
    methods (Access = public)
        function obj = SIM_Reconstructor(arg)
            %% 1st step: set default parameters in m_mParams
            obj.m_myParams = SetDefaultParams(obj.m_myParams);
            
            %% 2nd step: define all the commandline and config file options
            obj = setupProgramOptions(obj, arg);
            
            % If a config file is give, overwrite m_myParams from the config file
            if ~isempty(obj.m_config_file)
                try
                    obj.m_myParams = load(obj.m_config_file);
                catch
                    fprintf('can not open config file: %s\n', m_config_file);
                    fprintf('proceed without it\n');
                end
            end
            
            %% 3rd step: fill in m_myParams fields that have not been set yet
            obj = setParams(obj);
            fprintf('nphases=%d, ndirs=%d\n', obj.m_myParams.nphases, obj.m_myParams.ndirs);
            if obj.m_myParams.bNoKz0
                fprintf('Do not use kz=0 plane of the 0th order in the final assembly?  Yes, do not use.\n');
            else
                fprintf('Do not use kz=0 plane of the 0th order in the final assembly?  No, kz=0 plane is used.\n');
            end
            
            fprintf('DC singularity suppresion factor = %d\n', obj.m_myParams.dampenFactor);
            
            % TIFF mode:
            obj.m_all_matching_files = gatherMatchingFiles(obj.m_myParams.ifiles, obj.m_myParams.ofiles);
            % In TIFF mode, m_myParams.ifiles refers to the name of the folder raw data resides in;
            % and m_myParams.ofiles refers to a pattern in all the raw data file names.
            obj.m_imgParams.ntimes = size(obj.m_all_matching_files, 2);
            
            % save('config_file.mat', '-struct', 'obj.m_myParams');
            
            %% 4rd step: open m_myParams.ifiles; m_myParams.otffiles and check whether they exist
            obj = openFiles(obj);
            obj.m_zoffset = 0;
            
        end
        
        %% Load raw data from a file and do bleach correction
        function obj = loadAndRescaleImage(obj, timeIdx, waveIdx)
            %{
            * Calls private method loadImageData()
            * As is in other occasions, 'waveIdx' is to indicate color channel but mostly unused.
            * In the future, may add multi-color capability
            %}
            
            obj = loadImageData(obj, timeIdx, waveIdx);  % waveIdx is always 0 in most applications (single-color capability)
            
            if (obj.m_myParams.do_rescale)
                obj.m_reconData = rescaleDriver(timeIdx, waveIdx, obj.m_zoffset, obj.m_myParams, obj.m_imgParams, obj.m_reconData);
            end
            
            fprintf("Saving wide-field image for time point: #%d\n", timeIdx);
            saveWidefield(timeIdx, waveIdx, obj.m_myParams, obj.m_imgParams, obj.m_reconData.savedBands);
            fprintf("wide-field image saved\n\n");
            
            if (obj.m_myParams.bRLonOutput)
                fprintf("Saving wide-field Richardson-Lucy deconvolution image for time point: #%d\n", timeIdx);
                saveRLWidefield(timeIdx, waveIdx, obj.otf, obj.m_myParams, obj.m_imgParams, obj.m_reconData.savedBands);
                fprintf("Richardson-Lucy deconvolution wide-field image saved\n\n");
            end
            
            if (obj.m_myParams.bRLonInput)
                fprintf("2D Richardson-Lucy deconvolution on input images for time point: #%d\n", timeIdx);
                obj.m_reconData.savedBands = RL_inputdriver(obj.otf, obj.m_zoffset, obj.m_myParams, obj.m_imgParams, obj.m_reconData.savedBands);
                
                if (obj.m_myParams.do_rescale)
                    obj.m_reconData = rescaleDriver(timeIdx, waveIdx, obj.m_zoffset, obj.m_myParams, obj.m_imgParams, obj.m_reconData);
                end
                
                fprintf("Saving Richardson-Lucy deconvolution input images for time point: #%d\n", timeIdx);
                saveRLinput(timeIdx, waveIdx, obj.m_myParams, obj.m_imgParams, obj.m_reconData.savedBands);
                fprintf("Richardson-Lucy deconvolution input images saved\n\n");
            end
        end
        
        %% Return how many time points
        function ntimes = getNTimes(obj)
            ntimes = obj.m_imgParams.ntimes;
        end
        
        %% Set the current time index
        function obj = setCurTimeIdx(obj, it)
            obj.m_imgParams.curTimeIdx = it;
        end
        
        %% Process one SIM volume (i.e., for the current time point timeIdx)
        function obj = processOneVolume(obj, it)
            %{
            1. Re-initialize all key device buffers under m_reconData;
            2. Fine-tune k0 vectors, modulation amplitudes and phases for all directions and orders
            3. For all directions, pre-filter separated bands, inverse FFT, and assemble the bands
            %}
            
            obj.m_reconData.overlap0 = complex(zeros([obj.m_imgParams.ny, obj.m_imgParams.nx, obj.m_imgParams.nz], 'single'), zeros([obj.m_imgParams.ny, obj.m_imgParams.nx, obj.m_imgParams.nz], 'single'));
            obj.m_reconData.overlap1 = complex(zeros([obj.m_imgParams.ny, obj.m_imgParams.nx, obj.m_imgParams.nz], 'single'), zeros([obj.m_imgParams.ny, obj.m_imgParams.nx, obj.m_imgParams.nz], 'single'));
            
            [obj.m_myParams, obj.m_driftParams, obj.m_reconData] = findModulationVectorsAndPhasesForAllDirections(obj.m_zoffset, obj.otf, obj.m_myParams, obj.m_imgParams, obj.m_driftParams, obj.m_reconData);
            
            [obj.m_myParams, obj.m_imgParams] = calculateReconstructionParameters(obj.m_reconData.k0, obj.m_myParams, obj.m_imgParams);
            
            % Save memory space
            obj.m_reconData.overlap0 = [];
            obj.m_reconData.overlap1 = [];
           
            % RL part
            if (obj.m_myParams.bRLonOutput)
                obj.m_reconData.RLbuffer = zeros([obj.m_myParams.zoomfact * obj.m_imgParams.ny, obj.m_myParams.zoomfact * obj.m_imgParams.nx, obj.m_myParams.z_zoom * obj.m_imgParams.nz0], 'single');
                obj.m_reconData.SIM_OTF = zeros([obj.m_myParams.zoomfact * obj.m_imgParams.ny, obj.m_myParams.zoomfact * obj.m_imgParams.nx, obj.m_myParams.z_zoom * obj.m_imgParams.nz0], 'single');
                
                for direction = 1 : obj.m_myParams.ndirs
                    fprintf("RL deconvolution filter for angle %d\n\n", direction);
                    [shifted_result, otfSim] = filterbands_RL(direction, obj.m_reconData.savedBands(:, :, :, :, direction), obj.m_reconData.k0, obj.m_reconData.amp, obj.m_reconData.noiseVarFactors, obj.otf, obj.m_myParams, obj.m_imgParams);
                    
                    obj.m_reconData.RLbuffer = obj.m_reconData.RLbuffer + shifted_result;
                    obj.m_reconData.SIM_OTF = obj.m_reconData.SIM_OTF + otfSim;
                end
                
                ImageJ_formatted_TIFF.WriteTifStack(log(fftshift(abs(obj.m_reconData.RLbuffer))), 'RLbuffer.tif');
                obj.m_reconData.SIM_OTF = obj.m_reconData.SIM_OTF ./ obj.m_reconData.SIM_OTF(1, 1, 1);
                obj.m_reconData.RLbuffer = deconvolve(obj.m_reconData.RLbuffer, obj.m_reconData.SIM_OTF, obj.m_myParams.iter_output, true);
            end
            
            % Wiener part
            if (obj.m_myParams.botfBeforeShift)
                % loop all pattern directions
                for direction = 1 : obj.m_myParams.ndirs
                    fprintf("Wiener filtering for angle %d\n\n", direction);
                    filterbands_single_direction(it, direction, obj.m_reconData.savedBands(:, :, :, :, direction), ...
                        obj.m_reconData.k0, obj.m_reconData.amp, obj.m_reconData.noiseVarFactors, obj.otf, obj.m_myParams, obj.m_imgParams);
                end
                
            else
                fullResult = zeros([obj.m_myParams.zoomfact * obj.m_imgParams.ny, obj.m_myParams.zoomfact * obj.m_imgParams.nx, obj.m_myParams.z_zoom * obj.m_imgParams.nz0], 'single');
                % loop all pattern directions
                for direction = 1 : obj.m_myParams.ndirs
                    fprintf("Wiener filtering for angle %d\n\n", direction);
                    fullResult = fullResult + filterbands_single_direction_shifted(it, direction, obj.m_reconData.savedBands(:, :, :, :, direction), ...
                        obj.m_reconData.k0, obj.m_reconData.amp, obj.m_reconData.noiseVarFactors, obj.otf, obj.m_myParams, obj.m_imgParams);
                end
                % After loop all pattern directions, 'fullResult' now holds the image
            end
            
            
            obj.m_reconData.outbuffer = zeros([obj.m_myParams.zoomfact * obj.m_imgParams.ny, obj.m_myParams.zoomfact * obj.m_imgParams.nx, obj.m_myParams.z_zoom * obj.m_imgParams.nz0], 'single');
            if (obj.m_myParams.botfBeforeShift)
                % We will use Lin's method to do Wiener filtering                
                obj.m_reconData.bigbuffer = complex(zeros([obj.m_myParams.zoomfact * obj.m_imgParams.ny, obj.m_myParams.zoomfact * obj.m_imgParams.nx, obj.m_myParams.z_zoom * obj.m_imgParams.nz0], 'single'), ...
                    zeros([obj.m_myParams.zoomfact * obj.m_imgParams.ny, obj.m_myParams.zoomfact * obj.m_imgParams.nx, obj.m_myParams.z_zoom * obj.m_imgParams.nz0], 'single'));

                for direction = 1 : obj.m_myParams.ndirs
                    fprintf("Wiener filtering direction: #%d\n\n", direction);
                    obj.m_reconData.savedBands(:, :, :, :, direction) = filterbands(direction, obj.m_reconData.savedBands(:, :, :, :, direction), ...
                        obj.m_reconData.k0, obj.m_reconData.amp, obj.m_reconData.noiseVarFactors, obj.otf, obj.m_myParams, obj.m_imgParams);
                    
                    obj.m_reconData.outbuffer = assemblerealspacebands(direction, obj.m_reconData.outbuffer, obj.m_reconData.bigbuffer, ...
                        obj.m_reconData.savedBands(:, :, :, :, direction), obj.m_reconData.k0, obj.m_myParams, obj.m_imgParams);                    
                end
            else
                fprintf("Wiener filtering for all directions\n\n");
                % We will use 'fullResult' to do Wiener filtering
                obj.m_reconData.outbuffer = filterbands_shifted(obj.m_reconData.outbuffer, fullResult, obj.m_reconData.k0, obj.m_reconData.amp, obj.m_reconData.noiseVarFactors, obj.otf, obj.m_myParams, obj.m_imgParams);
            end
        end
        
        %% Write obj.m_reconData.outbuffer into tiff file
        function obj = writeResult(obj, it, iw)
            fprintf("Saving reconstructed results\n");
            
            if (obj.m_myParams.bRLonOutput)
                % Write RL SIM-OTF data image
                SIM_OTF = log(abs(fftshift(obj.m_reconData.SIM_OTF)));
                fileID = fopen(strcat('RL_SIM_OTF_' , num2str(it), '.tif'), 'w+');
                [Height, Width, Depth] = size(SIM_OTF);
                minval = min(SIM_OTF, [], 'all');
                maxval = max(SIM_OTF, [], 'all');
                header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(SIM_OTF), ...
                    1, Depth, 1, obj.m_imgParams.dz / obj.m_myParams.z_zoom, minval, maxval, obj.m_imgParams.dy / obj.m_myParams.zoomfact);
                SIM_OTF = permute(SIM_OTF, [2 1 3]);
                SIM_OTF = reshape(SIM_OTF, [1, Height * Width * Depth]);
                
                [~, ~, system_endian] = computer;
                if system_endian ~= header.endian
                    SIM_OTF = swapbytes(SIM_OTF);
                end
                
                SIM_OTF = typecast(SIM_OTF, 'uint8');
                fwrite(fileID, SIM_OTF, 'uint8');
                fclose(fileID);
                
                % Write RL deconvolution filtered image
                fileID = fopen(strcat('RL_' , num2str(it), '.tif'), 'w+');
                [Height, Width, Depth] = size(obj.m_reconData.RLbuffer);
                minval = min(obj.m_reconData.RLbuffer, [], 'all');
                maxval = max(obj.m_reconData.RLbuffer, [], 'all');
                header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(obj.m_reconData.RLbuffer), ...
                    1, Depth, 1, obj.m_imgParams.dz / obj.m_myParams.z_zoom, minval, maxval, obj.m_imgParams.dy / obj.m_myParams.zoomfact);
                RLbuffer = permute(obj.m_reconData.RLbuffer, [2 1 3]);
                RLbuffer = reshape(RLbuffer, [1, Height * Width * Depth]);
                
                if system_endian ~= header.endian
                    RLbuffer = swapbytes(RLbuffer);
                end
                
                RLbuffer = typecast(RLbuffer, 'uint8');
                fwrite(fileID, RLbuffer, 'uint8');
                fclose(fileID);
                
                % Wrtie RL deconvolution filtered fft data image
                RLbuffer_fft = log(abs(fftshift(fftn(obj.m_reconData.RLbuffer))));
                fileID = fopen(strcat('RL_fft_' , num2str(it), '.tif'), 'w+');
                [Height, Width, Depth] = size(RLbuffer_fft);
                minval = min(RLbuffer_fft, [], 'all');
                maxval = max(RLbuffer_fft, [], 'all');
                header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(RLbuffer_fft), ...
                    1, Depth, 1, obj.m_imgParams.dz / obj.m_myParams.z_zoom, minval, maxval, obj.m_imgParams.dy / obj.m_myParams.zoomfact);
                RLbuffer_fft = permute(RLbuffer_fft, [2 1 3]);
                RLbuffer_fft = reshape(RLbuffer_fft, [1, Height * Width * Depth]);
                
                if system_endian ~= header.endian
                    RLbuffer_fft = swapbytes(RLbuffer_fft);
                end
                
                RLbuffer_fft = typecast(RLbuffer_fft, 'uint8');
                fwrite(fileID, RLbuffer_fft, 'uint8');
                fclose(fileID);                
            end            
            

            if (obj.m_myParams.botfBeforeShift)
                % Write Wiener-filtered fft data image
                if (obj.m_myParams.bsavefft)
                    outbuffer_fft = obj.m_reconData.outbuffer;
                    if ~(obj.m_myParams.bfftkeepnegative)
                        outbuffer_fft(outbuffer_fft < 0) = 0;
                    end
                    outbuffer_fft = log(abs(fftshift(fftn(outbuffer_fft))));
                    fileID = fopen(strcat('Wiener_fft_' , num2str(it), '.tif'), 'w+');
                    [Height, Width, Depth] = size(outbuffer_fft);
                    minval = min(outbuffer_fft, [], 'all');
                    maxval = max(outbuffer_fft, [], 'all');
                    header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(outbuffer_fft), ...
                        1, Depth, 1, obj.m_imgParams.dz / obj.m_myParams.z_zoom, minval, maxval, obj.m_imgParams.dy / obj.m_myParams.zoomfact);
                    outbuffer_fft = permute(outbuffer_fft, [2 1 3]);
                    outbuffer_fft = reshape(outbuffer_fft, [1, Height * Width * Depth]);
                    
                    [~, ~, system_endian] = computer;
                    if system_endian ~= header.endian
                        outbuffer_fft = swapbytes(outbuffer_fft);
                    end
                    
                    outbuffer_fft = typecast(outbuffer_fft, 'uint8');
                    fwrite(fileID, outbuffer_fft, 'uint8');
                    fclose(fileID);
                end
                
                % Write Wiener-filtered image
                if ~(obj.m_myParams.bkeepnegative)
                    obj.m_reconData.outbuffer(obj.m_reconData.outbuffer < 0) = 0;
                end
                
                fileID = fopen(strcat(obj.m_myParams.ifiles, '\Wiener_' , num2str(it), '.tif'), 'w+');
                [Height, Width, Depth] = size(obj.m_reconData.outbuffer);
                minval = min(obj.m_reconData.outbuffer, [], 'all');
                maxval = max(obj.m_reconData.outbuffer, [], 'all');
                header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(obj.m_reconData.outbuffer), ...
                    1, Depth, 1, obj.m_imgParams.dz / obj.m_myParams.z_zoom, minval, maxval, obj.m_imgParams.dy / obj.m_myParams.zoomfact);
                outbuffer = permute(obj.m_reconData.outbuffer, [2 1 3]);
                outbuffer = reshape(outbuffer, [1, Height * Width * Depth]);
                
                [~, ~, system_endian] = computer;
                if system_endian ~= header.endian
                    outbuffer = swapbytes(outbuffer);
                end
                
                outbuffer = typecast(outbuffer, 'uint8');
                fwrite(fileID, outbuffer, 'uint8');
                fclose(fileID);
                
            else
                outbuffer = real(ifftn(obj.m_reconData.outbuffer));
                
                % Write Wiener-filtered fft data image
                if (obj.m_myParams.bsavefft)
                    outbuffer_fft = outbuffer;
                    if ~(obj.m_myParams.bfftkeepnegative)
                        outbuffer_fft(outbuffer_fft < 0) = 0;
                    end
                    outbuffer_fft = log(abs(fftshift(fftn(outbuffer_fft))));
                    
                    fileID = fopen(strcat('Wiener_fft_' , num2str(it), '.tif'), 'w+');
                    [Height, Width, Depth] = size(outbuffer_fft);
                    minval = min(outbuffer_fft, [], 'all');
                    maxval = max(outbuffer_fft, [], 'all');
                    header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(outbuffer_fft), ...
                        1, Depth, 1, obj.m_imgParams.dz / obj.m_myParams.z_zoom, minval, maxval, obj.m_imgParams.dy / obj.m_myParams.zoomfact);
                    outbuffer_fft = permute(outbuffer_fft, [2 1 3]);
                    outbuffer_fft = reshape(outbuffer_fft, [1, Height * Width * Depth]);
                    
                    [~, ~, system_endian] = computer;
                    if system_endian ~= header.endian
                        outbuffer_fft = swapbytes(outbuffer_fft);
                    end
                    
                    outbuffer_fft = typecast(outbuffer_fft, 'uint8');
                    fwrite(fileID, outbuffer_fft, 'uint8');
                    fclose(fileID);
                end
                
                % Write Wiener-filtered image
                if ~(obj.m_myParams.bkeepnegative)
                    outbuffer(outbuffer < 0) = 0;
                end
                
                fileID = fopen(strcat(obj.m_myParams.ifiles, '\Wiener_' , num2str(it), '.tif'), 'w+');
                [Height, Width, Depth] = size(outbuffer);
                minval = min(outbuffer, [], 'all');
                maxval = max(outbuffer, [], 'all');
                header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(outbuffer), ...
                    1, Depth, 1, obj.m_imgParams.dz / obj.m_myParams.z_zoom, minval, maxval, obj.m_imgParams.dy / obj.m_myParams.zoomfact);
                outbuffer = permute(outbuffer, [2 1 3]);
                outbuffer = reshape(outbuffer, [1, Height * Width * Depth]);
                
                [~, ~, system_endian] = computer;
                if system_endian ~= header.endian
                    outbuffer = swapbytes(outbuffer);
                end
                
                outbuffer = typecast(outbuffer, 'uint8');
                fwrite(fileID, outbuffer, 'uint8');
                fclose(fileID);
            end
            
            fprintf("Time point %d, wave %d done\n", it, iw);
        end
    end
    
    methods (Access = private)
        %% Define all the commandline and config file options
        function obj = setupProgramOptions(obj, arg)
            if ~isempty(arg.ifiles)
                obj.m_myParams.ifiles = arg.ifiles;
            else
                error("Input file name must be give.\n");
            end
            
            if ~isempty(arg.ofiles)
                obj.m_myParams.ofiles = arg.ofiles;
            else
                error("Output file name must be give.\n");
            end
            
            if ~isempty(arg.otffiles)
                obj.m_myParams.otffiles = arg.otffiles;
            else
                error("OTF file name must be give.\n");
            end
            
            if ~isempty(arg.usecorr)
                obj.m_myParams.corrfiles = arg.usecorr;
            end
            
            if ~isempty(arg.ndirs)
                obj.m_myParams.ndirs = arg.ndirs;
            else
                obj.m_myParams.ndirs = 3;
            end
            
            if ~isempty(arg.nphases)
                obj.m_myParams.nphases = arg.nphases;
            else
                obj.m_myParams.nphases = 5;
            end
            
            if ~isempty(arg.nordersout)
                obj.m_myParams.norders_output = arg.nordersout;
            else
                obj.m_myParams.norders_output = 0;
            end
            
            if ~isempty(arg.angle0)
                obj.m_myParams.k0startangle = arg.angle0;
            else
%                 obj.m_myParams.k0startangle = -0.008;
                obj.m_myParams.k0startangle = [];
            end
            
            if ~isempty(arg.ls)
                obj.m_myParams.linespacing = arg.ls;
            else
%                 obj.m_myParams.linespacing = 0.198802;
                obj.m_myParams.linespacing = [];
            end
            
            if ~isempty(arg.na)
                obj.m_myParams.na = arg.na;
            else
                obj.m_myParams.na = 1.35;
            end
            
            if ~isempty(arg.nimm)
                obj.m_myParams.nimm = arg.nimm;
            else
                obj.m_myParams.nimm = 1.406;
            end
            
            if ~isempty(arg.nmed)
                obj.m_myParams.nmed = arg.nmed;
            else
                obj.m_myParams.nmed = 1.333;
            end
            
            if ~isempty(arg.iter_input)
                obj.m_myParams.iter_input = arg.iter_input;
            else
                obj.m_myParams.iter_input = 1;
            end
            
            if ~isempty(arg.iter_output)
                obj.m_myParams.iter_output = arg.iter_output;
            else
                obj.m_myParams.iter_output = 5;
            end
            
            if ~isempty(arg.zoomfact)
                obj.m_myParams.zoomfact = arg.zoomfact;
            else
                obj.m_myParams.zoomfact = 2;
            end
            
            if ~isempty(arg.explodefact)
                obj.m_myParams.explodefact = arg.explodefact;
            else
                obj.m_myParams.explodefact = 1.0;
            end
            
            if ~isempty(arg.zzoom)
                obj.m_myParams.z_zoom = arg.zzoom;
            else
                obj.m_myParams.z_zoom = 1;
            end
            
            if isfield(arg, 'nofilteroverlaps')
                obj.m_myParams.bFilteroverlaps = false;
            end
            
            if ~isempty(arg.background)
                obj.m_myParams.constbkgd = arg.background;
            else
                obj.m_myParams.constbkgd = 0.0;
            end
            
            if ~isempty(arg.wiener)
                obj.m_myParams.wiener = arg.wiener;
            else
                obj.m_myParams.wiener = 0.01;
            end
            
            if ~isempty(arg.linear_wiener)
                obj.m_myParams.linearWiener = arg.linear_wiener;
            else
                obj.m_myParams.linearWiener = 0;
            end
            
            if ~isempty(arg.forcemodamp)
                obj.m_myParams.forceamp = arg.forcemodamp;
            end
            
            if ~isempty(arg.k0angles)
                obj.m_myParams.k0angles = arg.k0angles;
            end
            
            if isfield(arg, 'useRLonInput')
                obj.m_myParams.bRLonInput = true;
            end
            
            if isfield(arg, 'useRLonOutput')
                obj.m_myParams.bRLonOutput = true;
            end
            
            if isfield(arg, 'otfRA')
                obj.m_myParams.bRadAvgOTF = true;
            end
            
            if isfield(arg, 'bOneOTFperAngle')
                obj.m_myParams.bRadAvgOTF = true;
            end
            
            if isfield(arg, 'fastSI')
                obj.m_myParams.bFastSIM = true;
            end
            
            if isfield(arg, 'searchforvector')
                obj.m_myParams.bSearchforvector = true;
            end
                
            if isfield(arg, 'k0searchAll')
                obj.m_myParams.bUseTime0k0 = false;
            end
            
            if isfield(arg, 'equalizez')
                obj.m_myParams.equalizez = true;
            end
            
            if isfield(arg, 'equalizet')
                obj.m_myParams.equalizet = true;
            end
            
            if isfield(arg, 'dampenOrder0')
                obj.m_myParams.bDampenOrder0 = true;
            end
            
            if isfield(arg, 'nosuppress')
                obj.m_myParams.bSuppress_singularities = false;
            end
            
            if isfield(arg, 'dampenFactor')
                obj.m_myParams.dampenFactor = arg.dampenFactor;
            end
            
            if isfield(arg, 'nokz0')
                obj.m_myParams.bNoKz0 = true;
            end
            
            if isfield(arg, 'applyOtfBeforeShift')
                obj.m_myParams.botfBeforeShift = true;
            else
                obj.m_myParams.botfBeforeShift = false;
            end
            
            if ~isempty(arg.gammaApo)
                obj.m_myParams.apoGamma = arg.gammaApo;
            else
                obj.m_myParams.apoGamma = 1.0;
            end
            
            if ~isempty(arg.saveprefiltered)
                obj.m_myParams.fileSeparated = arg.saveprefiltered;
            end
            
            if ~isempty(arg.savealignedraw)
                obj.m_myParams.fileRawAligned = arg.savealignedraw;
            end
            
            if ~isempty(arg.saveoverlaps)
                obj.m_myParams.fileOverlaps = arg.saveoverlaps;
            end
            
            if ~isempty(arg.config)
                obj.m_config_file = arg.config;
            else
                obj.m_config_file = [];
            end
            
            if isfield(arg, 'twolenses')
                obj.m_myParams.bTwolens = true;
            end
            
            if ~isempty(arg.lambda)
                obj.m_myParams.lambda = arg.lambda;
            else
                obj.m_myParams.lambda = 488;
            end            
            
            if ~isempty(arg.wavelength)
                obj.m_imgParams.wave(1) = arg.wavelength;
            else
                obj.m_imgParams.wave(1) = 525;
            end
            
            if isfield(arg, 'keepnegativevalues')
                obj.m_myParams.bkeepnegative = true;
            else
                obj.m_myParams.bkeepnegative = false;
            end
            
            if isfield(arg, 'savefft')
                obj.m_myParams.bsavefft = true;
            else
                obj.m_myParams.bsavefft = false;
            end
            
            if isfield(arg, 'fftkeepnegativevalues')
                obj.m_myParams.bfftkeepnegative = true;
            else
                obj.m_myParams.bfftkeepnegative = false;
            end
        end
        
        %% Fill in m_myParams fields that have not been set yet
        function obj = setParams(obj)
            if ~isempty(obj.m_myParams.corrfiles)
                obj.m_myParams.bUsecorr = true;
            end
            
            if obj.m_myParams.wiener > 0
                obj.m_myParams.bUseEstimatedWiener = false;
            end
            if obj.m_myParams.linearWiener > 0
                fprintf('wiener=%.3f (+%.5f)\n', obj.m_myParams.wiener, obj.m_myParams.linearWiener);
            else
                fprintf('wiener=%.3f\n', obj.m_myParams.wiener);
            end
            
            
            if ~isempty(obj.m_myParams.apoGamma)
                obj.m_myParams.apodizeoutput = 2;
                fprintf('gamma=%.1f\n', obj.m_myParams.apoGamma);
            end
            
            if ~isempty(obj.m_myParams.fileSeparated)
                obj.m_myParams.bSaveSeparated = true;
            end
            
            if ~isempty(obj.m_myParams.fileRawAligned)
                obj.m_myParams.bSaveAlignedRaw = true;
            end
            
            if ~isempty(obj.m_myParams.fileOverlaps)
                obj.m_myParams.bSaveOverlaps = true;
            end
            
            if ~isempty(obj.m_myParams.k0angles)
                fprintf('k0angles=%.3f %.3f %.3f\n', obj.m_myParams.k0angles);
            end
        end
        
        %% Open m_myParams.ifiles; m_myParams.otffiles and check whether they exist
        function obj = openFiles(obj)
            % TIFF mode:
            % In TIFF mode, input files are not opened till loadAndRescaleImage()
            if ~size(obj.m_all_matching_files, 1)
                error('Input file not found\n');
            end
            
            % In TIFF mode, output files are not created until writeResult() is called
            [otfstream_no, errmsg] = fopen(obj.m_myParams.otffiles, 'r');
            if not(isempty(errmsg))
                warning(errmsg);
                error('OTF file not found\n');
            end
            fclose(otfstream_no);
            
        end
        
        %% Load raw data from a file and do bleach correction (what "rescale" refers to)
        function obj = loadImageData(obj, it, iw)
            [~, ~, system_endian] = computer;
            
            % TIFF mode:
            raw_tiff = obj.m_all_matching_files{it + 1};
            % set up m_myParams, m_imgParams, and m_reconData based on the first input TIFF
            if it == 0
                obj = setup(obj, raw_tiff);
            end
            
            if (obj.m_imgParams.nz0 < obj.m_imgParams.nz)
                nz = obj.m_imgParams.nz0;
            else
                nz = obj.m_imgParams.nz;
            end
            
            istream_no = fopen(raw_tiff, 'r');
            ifile_header = ImageJ_formatted_TIFF.parse_tif(istream_no, 0);
            
            % Calculate data_type from header.BitsPerSample
            data_type = [ifile_header.SampleFormat, ifile_header.BitsPerSample];
            if isequal(data_type, [1, 8])
                dtype = 'uint8';
            elseif isequal(data_type, [1, 16])
                dtype = 'uint16';
            elseif isequal(data_type, [1, 32])
                dtype = 'uint32';
            elseif isequal(data_type, [1, 64])
                dtype = 'uint64';
            elseif isequal(data_type, [2, 8])
                dtype = 'int8';
            elseif isequal(data_type, [2, 16])
                dtype = 'int16';
            elseif isequal(data_type, [2, 32])
                dtype = 'int32';
            elseif isequal(data_type, [2, 64])
                dtype = 'int64';
            elseif isequal(data_type, [3, 32])
                dtype = 'single';
            elseif isequal(data_type, [3, 64])
                dtype = 'double';
            else
                error("ReadTifStack does not support SampleFormat:%d with BitsPerSample: %d", ifile_header.SampleFormat, ifile_header.BitsPerSample);
            end
            
            % Initialize reading buffer and parameters
            PixelNum = obj.m_imgParams.nx_input * obj.m_imgParams.ny_input;
            ByteCounts = fix(PixelNum * ifile_header.BitsPerSample / 8);
            buffer = zeros([1, PixelNum], dtype);
            image_input = zeros([obj.m_imgParams.ny_input, obj.m_imgParams.nx_input, obj.m_imgParams.nz0, obj.m_myParams.nphases, obj.m_myParams.ndirs], 'single');
            fseek(istream_no, ifile_header.StripOffsets, 'bof');
            
            if (obj.m_myParams.bFastSIM)
                % data organized into (nz, ndirs, nphases)
                for z = 1 : nz
                    for direction = 1 : obj.m_myParams.ndirs
                        for phase = 1 : obj.m_myParams.nphases
                            buffer(:) =  typecast(transpose(fread(istream_no, ByteCounts, 'uint8=>uint8')), dtype);
                            if system_endian ~= ifile_header.endian
                                buffer = swapbytes(buffer);
                            end
                            image_input(:, :, z + obj.m_zoffset, phase, direction) = single(permute(reshape(buffer, [obj.m_imgParams.nx_input, obj.m_imgParams.ny_input]), [2 1]));
                        end
                    end
                end
            else
                % data organized into (ndirs, nz, nphases)
                for direction = 1 : obj.m_myParams.ndirs
                    for z = 1 : nz
                        for phase = 1 : obj.m_myParams.nphases
                            buffer(:) =  typecast(transpose(fread(istream_no, ByteCounts, 'uint8=>uint8')), dtype);
                            if system_endian ~= ifile_header.endian
                                buffer = swapbytes(buffer);
                            end
                            image_input(:, :, z + obj.m_zoffset, phase, direction) = single(permute(reshape(buffer, [obj.m_imgParams.nx_input, obj.m_imgParams.ny_input]), [2 1]));
                        end
                    end
                end
            end
            
            fclose(istream_no);
            
            if (obj.m_myParams.bPadBlur)
                for direction = 1 : obj.m_myParams.ndirs
                    for phase = 1 : obj.m_myParams.nphases
                        dkr = 1.0 / (obj.m_imgParams.ny * obj.m_imgParams.dy);
                        krscale = dkr / obj.otf.dkrotf;
                        otf_2D = obj.otf.getOtfMatrix2D(0, 1, obj.m_imgParams.nx, obj.m_imgParams.ny, 0, 0, krscale);
                        blur_top = abs(ifftn(fftn(image_input(:, :, 1 + obj.m_zoffset, phase, direction)) .* otf_2D));
                        blur_bottom = abs(ifftn(fftn(image_input(:, :, nz + obj.m_zoffset, phase, direction)) .* otf_2D));
                        for z = 1 : obj.m_zoffset
                            image_input(:, :, z, phase, direction) = blur_top;
                            image_input(:, :, z + nz + obj.m_zoffset, phase, direction) = blur_bottom;
                        end
                    end
                end                
            end
            
            % Subtracts a constant value (background) from each pixel
            image_input = image_input - obj.m_myParams.constbkgd;
            image_input(image_input < 0) = 0;
            
            % Apodize (or edge softening) every 2D slice:
            image_input = apodizationDriver(obj.m_zoffset, obj.m_myParams, obj.m_imgParams, image_input);

            % Clear obj.m_reconData.savedBands
            obj.m_reconData.savedBands(:) = 0;
            % Put edge-softened, background-subtracted raw images to savedBands
            obj.m_reconData.savedBands(1 : obj.m_imgParams.ny_input, 1 : obj.m_imgParams.nx_input, :, :, :) = image_input;

            % inscale is used to normalize all image layers so that their integral is 1
%             obj.m_reconData.savedBands = obj.m_reconData.savedBands * obj.m_imgParams.inscale;
            
        end
        
        %% Read header, set up OTF and separation matrix, allocate device buffers, initialize parameters like k0guess etc.
        function obj = setup(obj, inTIFF)
            istream_no = fopen(inTIFF, 'r');
            ifile_header = ImageJ_formatted_TIFF.parse_tif(istream_no, 0);
            fclose(istream_no);
            
            obj.m_imgParams.nx_input = ifile_header.ImageWidth;     % x pixel number of input image
            obj.m_imgParams.ny_input = ifile_header.ImageLength;    % y pixel number of input image
            
            % header.nz is defined as the total number of sections = nz * nwaves * ntimes * (ndirs * nphases) in our case
            ImageDescription = convertStringsToChars(ifile_header.ImageDescription);
            % Calculate nz from header.ImageDescription
            k1 = strfind(ImageDescription, 'images=');
            ImageDescription_crop = ImageDescription(k1:end);
            k2 = strfind(ImageDescription_crop, newline);
            if ~isempty(k1) && ~isempty(k2)
                depth = str2double(ImageDescription(k1(1) + 7: k1(1) + k2(1) - 2));
            else
                warning("Did not find 'images=' in ImageDescription. Try to calculate it from filesize.");
                FileAttributes = dir(inTIFF);
                depth = floor(FileAttributes.bytes / (ifile_header.ImageWidth * ifile_header.ImageLength * (ifile_header.BitsPerSample / 8)));
            end
            obj.m_imgParams.nz = depth;
            obj.m_imgParams.nwaves = 1;   % Because we save each emission wavelength as individual image, imgParams.nwaves is always 1
            obj.m_imgParams.nz = fix(obj.m_imgParams.nz / (obj.m_myParams.nphases * obj.m_myParams.ndirs));  % z pixel number of input image
            
            % Calculate spacing from header.
            k1 = strfind(ImageDescription, 'spacing=');
            ImageDescription_crop = ImageDescription(k1:end);
            k2 = strfind(ImageDescription_crop, newline);
            if ~isempty(k1) && ~isempty(k2)
                spacing = str2double(ImageDescription(k1(1) + 8: k1(1) + k2(1) - 2));
            else
                warning("Did not find 'spacing=' in ImageDescription. Suppose it is 0.04 microns.");
                spacing = 0.04;
            end
            % dy: lateral pixel size; dz: axial pixel size; both in microns
            obj.m_imgParams.dy = ifile_header.resolution;
            obj.m_imgParams.dz = spacing;
            % obj.m_imgParams.wave(1) has been imported from setupProgramOptions()
            
            if obj.m_myParams.nzPadTo
                obj.m_imgParams.nz0 = obj.m_myParams.nzPadTo;
            else
                obj.m_imgParams.nz0 = obj.m_imgParams.nz;
            end
            fprintf("Input images: nx=%d, ny=%d, nz=%d, nz0=%d, nwaves=%d, ntimes=%d\n", obj.m_imgParams.nx_input, obj.m_imgParams.ny_input, obj.m_imgParams.nz, obj.m_imgParams.nz0, obj.m_imgParams.nwaves, obj.m_imgParams.ntimes);
            
            % to make sure the image is a square with multiples of 32 instead of rectangle
            if ((obj.m_imgParams.nx_input ~= obj.m_imgParams.ny_input) || (mod(obj.m_imgParams.nx_input, 32) ~= 0) || (mod(obj.m_imgParams.ny_input, 32) ~= 0))
                imgSize = 32 * ceil(max(obj.m_imgParams.nx_input, obj.m_imgParams.ny_input) / 32);
                % Re-define (m_imgParams.nx, m_imgParams.ny) if we have to resize the images...
                obj.m_imgParams.nx = imgSize;
                obj.m_imgParams.ny = imgSize;
                fprintf("Resizing input images to: nx=%d, ny=%d\n", imgSize, imgSize);
            else
                % Or just take the same width & height as input images
                obj.m_imgParams.nx = obj.m_imgParams.nx_input;
                obj.m_imgParams.ny = obj.m_imgParams.ny_input;
            end
            
            [obj.m_myParams, obj.m_imgParams, obj.m_reconData, obj.otf] = setup_part2(obj.m_myParams, obj.m_imgParams, obj.m_reconData);
            
            if (obj.m_myParams.nzPadTo)
                obj.m_zoffset = fix((obj.m_imgParams.nz0 - obj.m_imgParams.nz) / 2);
            end
            
            % Calculate other important image parameters, including dkr, dkz, krscale, kzscale
            
            % dkr: lateral frequency / pixel in reciprocal space;
            obj.m_imgParams.dkr = 1.0 / (obj.m_imgParams.ny * obj.m_imgParams.dy);
            
            % dkz: axial frequency / pixel in reciprocal space;
            if (obj.m_imgParams.dz > 0)
                obj.m_imgParams.dkz = 1.0 / (obj.m_imgParams.nz0 * obj.m_imgParams.dz);
            else
                obj.m_imgParams.dkz = obj.otf.dkzotf;
            end
            
            % ratio of radial direction pixel scales of data and otf
            obj.m_imgParams.krscale = obj.m_imgParams.dkr / obj.otf.dkrotf;
            
            % ratio of axial direction pixel scales of data and otf
            obj.m_imgParams.kzscale = obj.m_imgParams.dkz / obj.otf.dkzotf;
        end
        
    end
end





%% Set default parameters in m_mParams
function Params = SetDefaultParams(Params)
Params.k0startangle = -0.008;
Params.linespacing =  0.198802;         % default to Olympus 1.35N.A. 100x silicone oil
Params.na = 1.35;
Params.nimm = 1.406;
Params.nmed =1.333;
Params.ndirs = 3;
Params.nphases = 5;
Params.norders_output  = 0;
Params.phaseSteps = 0;
Params.bTwolens = false;
Params.bFastSIM = false;
Params.bRLonInput = false;
Params.iter_input = 1;
Params.bRLonOutput = false;
Params.iter_output = 5;
Params.zoomfact = 2;
Params.z_zoom = 1;
Params.nzPadTo = 0;
Params.bPadBlur = false;
Params.explodefact = 1.0;
Params.bFilteroverlaps = true;
Params.recalcarrays = 1;
Params.napodize = 10;
% Params.k0angles = [];
Params.bSearchforvector = false; % Pay attention here.
Params.bUseTime0k0 = true;
Params.apodizeoutput = 0;
Params.apoGamma = 1.0;
Params.bSuppress_singularities = true;
Params.dampenFactor = 1000;
Params.suppression_radius = 10;
Params.bDampenOrder0 = false;
Params.bFitallphases = true;
Params.do_rescale = 1;
Params.equalizez = false;
Params.equalizet = false;
Params.bNoKz0 = false;
Params.botfBeforeShift = false;
Params.wiener = 0.01;
Params.bUseEstimatedWiener = true;
Params.bRadAvgOTF = false;
Params.bOneOTFperAngle = false;
Params.bFixdrift = false;
Params.drift_filter_fact = 0.0;
Params.constbkgd = 0.0;
Params.bBgInExtHdr = false;
Params.bUsecorr = false;
Params.bMakemodel = false;
Params.bSaveSeparated = false;
Params.bSaveAlignedRaw = false;
Params.bSaveOverlaps = false;
Params.lambda = 488;
end

%% In TIFF mode, find all the tiff file with given filename pattern
function matchingfiles = gatherMatchingFiles(target_path, pattern)
matchingfiles = {};
allfiles = dir([target_path '\*.tif']);
for i = 1 : size(allfiles, 1)
    if contains(allfiles(i).name, pattern)
        matchingfiles{end + 1} = strcat(allfiles(i).folder, '\', allfiles(i).name); % Dynamic, cannot be preallocated.
    end
end
end

%% Allocate memory & initialize a bunch of reconstruction related parameters
function [params, imgParams, reconData, otf] = setup_part2(params, imgParams, reconData)
%{
1. OTF is loaded and transferred to RAM by calling getOTFs();
2. Separation (or un-mixing) matrix is allocated and populated;
3. Raw image and other buffers, including overlap0/1, bigbuffer, outbuffer,
   are allocated on ram by calling allocateImageBuffers();
4. Allocate/initialize a bunch of reconstruction related parameters , including:
 a. inscale
 b. k0
 c. k0_time0
 d. k0guess
 e. sum_dir0_phase0
 f. amp
%}

%% 1st step: OTF is loaded and transferred to RAM;
[params, otf] = getOTFs(params, imgParams);

% To allocate Separation Matrix and Noise Variance Amplification Factors
reconData = allocSepMatrixAndNoiseVarFactors(params, reconData);

%% 2nd step: generates the matrix that is to be used to separate the raw indata into the different bands of sample information
[reconData.sepMatrix, reconData.noiseVarFactors] = makematrix(params.nphases, params.norders, 1, 0, reconData.sepMatrix, reconData.noiseVarFactors);

%% 3rd step: Raw image and other buffers, including overlap0/1, bigbuffer, outbuffer, are allocated on RAM by calling allocateImageBuffers()
reconData = allocateImageBuffers(params, imgParams, reconData);

%% 4th step: Allocate/initialize a bunch of reconstruction related parameters, including: a.inscale; b.k0; c.k0_time0; d.k0guess; e.sum_dir0_phase0; f.amp

% % Here, I overwrite obj.m_myParams.napodize by 6% of rim width,
% % or you can always use obj.m_myParams.napodize = 10 pixels;
% if (params.napodize >= 0) && (params.napodize < fix(min(imgParams.nx_input, imgParams.ny_input) * 0.06))
%     params.napodize = fix(min(imgParams.nx_input, imgParams.ny_input) * 0.06);
% end

imgParams.inscale = 1.0 / (imgParams.nx * imgParams.ny * imgParams.nz0 * params.zoomfact * params.zoomfact * params.z_zoom * params.ndirs);

k0guess = struct('x',[], 'y', []);
reconData.k0 = repmat(k0guess, 1, params.ndirs);
reconData.k0_time0 = repmat(k0guess, 1, params.ndirs);
reconData.k0guess = repmat(k0guess, 1, params.ndirs);

if (~isempty(params.k0startangle) && ~isempty(params.linespacing))
    delta_angle = pi / params.ndirs;                    % delta_angle = pi / 3
    dkr = 1 / (imgParams.ny * imgParams.dy);            % dkr = 1 / (1024 * 0.039) = 0.0250 um-1
%     k0magguess = (1.0 / params.linespacing) / dkr;      % k0magguess = 1.0 / 0.198802 / 0.0250 = 200.8833
%     if (imgParams.nz > 1)
%         % 3D-SIM case
%         nordersIn = fix(params.nphases / 2) + 1;
%         k0magguess = k0magguess / (nordersIn - 1);      % k0magguess = 100.4416
%     elseif (imgParams.nz == 1)
%         % 2D-SIM case: Here we are using 5 phases to do 2D-SIM, so it is still:
%         nordersIn = fix(params.nphases / 2) + 1;
%         k0magguess = k0magguess / (nordersIn - 1);      % k0magguess = 100.4416
%     end
    for i = 0 : (params.ndirs - 1)
        if (size(params.k0angles, 2) < params.ndirs)
            if i == 0
                k0angleguess = params.k0startangle;
            elseif i == 1
                k0angleguess = params.k0startangle - delta_angle;
            elseif i == 2
                k0angleguess = params.k0startangle + delta_angle;
            end
        else
            k0angleguess = params.k0angles(i + 1);
        end

        k0magguess = (1.0 / params.linespacing(i + 1)) / dkr;
        if (imgParams.nz > 1)
            % 3D-SIM case
            nordersIn = fix(params.nphases / 2) + 1;
            k0magguess = k0magguess / (nordersIn - 1);      % k0magguess = 100.4416
        elseif (imgParams.nz == 1)
            % 2D-SIM case: Here we are using 5 phases to do 2D-SIM, so it is still:
            nordersIn = fix(params.nphases / 2) + 1;
            k0magguess = k0magguess / (nordersIn - 1);      % k0magguess = 100.4416
        end

        k0guess.x = k0magguess * cos(k0angleguess);
        k0guess.y = k0magguess * sin(k0angleguess);
        reconData.k0guess(i + 1) = k0guess;
    end
end

reconData.sum_dir0_phase0 = zeros([imgParams.nz * imgParams.nwaves, 1]);
reconData.amp = complex(zeros([params.norders, params.ndirs], 'single'), zeros([params.norders, params.ndirs], 'single'));

for i = 1 : params.ndirs
    reconData.amp(1, i) = complex(single(1.0), single(0.0));    % Amp is always 1 for the 0th order
end

end

%% 1st step: OTF is loaded and transferred to RAM by calling getOTFs()
function [params, otf] = getOTFs(params, imgParams)

params.norders = 0;
if params.norders_output ~= 0
    params.norders = params.norders_output;
else
    params.norders = fix(params.nphases / 2) + 1;
end

% Create an OTF object and transfer params
otf = OtfProvider(params, imgParams);

% Get the OTF header information from otf-file. So (nxotf, nyotf, nzotf, dkrotf, dkzotf) are defined in obj.otk
otf = otf.determine_otf_dimensions(imgParams);

% Transfer OTF data from .mat file to RAM
otf = otf.loadOTFs();

end

%% Allocate Separation Matrix and Noise Variance Amplification Factors
function reconData = allocSepMatrixAndNoiseVarFactors(params, reconData)
reconData.sepMatrix = zeros([params.norders * 2 - 1, params.nphases], 'single');
reconData.noiseVarFactors = ones([params.ndirs, params.norders], 'single');
end

%% 2nd step: Generates the matrix that is to be used to separate the raw indata into the different bands of sample information
function [sepMatrix, noiseVarFactors] = makematrix(nphases, norders, dire, arrPhases, sepMatrix, noiseVarFactors)
%{
Two cases:
	1. If arrPhases is NULL, we're certain that phases are equally spaced within 2Pi range;
	2. Else, user supplies a list of phases actually used (e.g, when there is drift or
	non-ideal SLM patterns are used), and noise variance amplification factors for each
	order and direction is calculated based on the non-ideal separation matrix
%}

%    phase = 0;         2pi / 5					2 * 2pi / 5				3 * 2pi / 5				4 * 2pi / 5
% sepMatrix = [1,			1,						1,						1,					1;                      order = 0
%			  1,	cos(1 * 1 * 2pi / 5),   cos(2 * 1 * 2pi / 5),   cos(3 * 1 * 2pi / 5),	cos(4 * 1 * 2pi / 5);		order = 1, bandre
%			  0,	sin(1 * 1 * 2pi / 5),   sin(2 * 1 * 2pi / 5),   sin(3 * 1 * 2pi / 5),	sin(4 * 1 * 2pi / 5);		order = 1, bandim
%			  1,	cos(1 * 2 * 2pi / 5),   cos(2 * 2 * 2pi / 5),   cos(3 * 2 * 2pi / 5),	cos(4 * 2 * 2pi / 5);		order = 2, bandre
%			  0,	sin(1 * 2 * 2pi / 5),   sin(2 * 2 * 2pi / 5),   sin(3 * 2 * 2pi / 5),	sin(4 * 2 * 2pi / 5);]		order = 2, bandim

fprintf("In makematrix.\n");

if arrPhases == 0
    phi = 2.0 * pi / nphases;
    for j = 0 : (nphases - 1)
        sepMatrix(1, j + 1) = 1.0;
        for order = 1 : (norders - 1)
            %{
            With this definition:   bandplus = bandre - i bandim;
                                    bandminus = bandre + i bandim
            has coefficient exp() with unit normalization
            %}
            sepMatrix(2 * order, j + 1) = cos(j * order * phi);
            sepMatrix(2 * order + 1, j + 1) = sin(j * order * phi);
        end
        
    end
else
    %{
    Here we differentiate between 2 options:
		/* 1. use direct matrix inversion;*/
		/* 2. use least-square solution.*/
    %}
    nCols = 2 * norder - 1;
    if (nCols < nphases)
        % use least-square
        forwardM = zeros([nCols, params.nphases], 'single');
        A_t_A = zeros([nCols, nCols], 'single');
        % First construct the forward matrix
        for i = 0 : (nphases - 1)
            forwardM(1, i + 1) = 1.0  / nphases;
            for order = 1 : (norders - 1)
                forwardM(2 * order, i + 1) = 2 * cos(order * arrPhases(i + 1)) / nphases;
                forwardM(2 * order + 1, i + 1) = 2 * sin(order * arrPhases(i + 1)) / nphases;
            end
        end
        
        % multiply transposed forward matrix with forward matrix
        A_t_A(:) = forwardM * transpose(forwardM);
        % Then invert the forward matrix to form the sep matrix
        %         A_t_A = inv(A_t_A);
        % multiply inverse A_t_A with transposed forward matrix
        sepMatrix = transpose(forwardM) / A_t_A;
        % transpose back to C-style matrix
        sepMatrix = transpose(sepMatrix);
    else
        % use direct inversion
        % First construct the forward matrix (in C-style storage convention )
        for i = 0 : (nphases - 1)
            sepMatrix(1, i + 1) = 1.0 / nphases;
            for order = 1 : (norders - 1)
                sepMatrix(2 * order, i + 1) = 2 * cos(order * arrPhases(i + 1)) / nphases;
                sepMatrix(2 * order + 1, i + 1) = 2 * sin(order * arrPhases(i + 1)) / nphases;
            end
        end
        
        % Then invert the forward matrix to form the sep matrix
        sepMatrix = inv(sepMatrix);
    end
    
    % Report noise factors
    for order = 0 : (norders - 1)
        noiseVarFactors(dire, order + 1) = 0.0;
    end
    
    for j = 0 : (nphases - 1)
        noiseVarFactors(dire, 1) = noiseVarFactors(dire, 1) + sepMatrix(1, j + 1) ^ 2;
        for order = 1 : (norders - 1)
            noiseVarFactors(dire, order) = noiseVarFactors(dire, order) + sepMatrix(2 * order, j + 1) ^ 2 + sepMatrix(2 * order + 1, j + 1) ^ 2;
        end
    end
    
    fprintf(" Separation noise factors: ");
    for order = 0 : (norders - 1)
        noiseVarFactors(dire, order + 1) = noiseVarFactors(dire, order + 1) / nphases;
        fprintf(" Order %d: %.3f,", order, sqrt(noiseVarFactors(dir, order + 1)));
    end
    fprintf('\n');
end

fprintf("Separation matrix:\n");
for j = 1 : (norders * 2 - 1)
    for i = 1 : nphases
        fprintf("%9.5f ", sepMatrix(j, i));
    end
    fprintf("\n");
end
fprintf("\n");
end

%% 3rd step: Allocate RAM for m_reconData.savedBands and set all to zero
function data = allocateImageBuffers(params, imgParams, data)
data.savedBands = zeros([imgParams.ny, imgParams.nx, imgParams.nz0, params.nphases, params.ndirs], 'single');
end

%% Apodize (or edge softening) every 2D slice
function data = apodizationDriver(zoffset, params, imgParams, data)

for direction = 1 : params.ndirs
    rawImages = data(:, :, :, :, direction);    
    for z = 1 : imgParams.nz
        for phase = 1 : params.nphases
            if (params.napodize >= 0)
                data(:, :, z + zoffset, phase, direction) = apodize(params.napodize, imgParams.nx_input, imgParams.ny_input, rawImages(:, :, z + zoffset, phase));
            elseif (params.napodize == -1)
                data(:, :, z + zoffset, phase, direction) = cosapodize(params.napodize, imgParams.nx_input, imgParams.ny_input, rawImages(:, :, z + zoffset, phase));
            end
        end
    end
end

end

%% to apodize (or edge softening) every 2D slice
function image = apodize(napodize, nx, ny, image)

% Lin's algorithm

% % Row-wise soften
% diff = (image(ny, :) - image(1, :)) / 2;    % Difference between last (512) row and the first (1) row;
% for i = 1: napodize
%     fact = 1 - sin(((i - 0.5) / napodize) * pi * 0.5);
%     image(i, :) = image(i, :) + diff .* fact;
%     image(ny + 1 - i, :) = image(ny + 1 - i, :) - diff .* fact;
% end
% 
% % Column-wise soften
% diff = (image(:, nx) - image(:, 1)) / 2;   % Difference between last (256) column and the first (1) column;
% for j = 1: napodize
%     fact = 1 - sin(((j - 0.5) / napodize) * pi * 0.5);
%     image(:, j) = image(:, j) + diff .* fact;
%     image(:, nx + 1 - j) = image(:, nx + 1 - j) - diff .* fact;
% end

% My algorithm # 1

% % Row-wise soften
% diff1 = 0 - image(1, :);    % Difference between 0 and the first (1) row;
% diff2 = 0 - image(ny, :);   % Difference between 0 and the last (512) row;
% for i = 1: napodize
%     fact = 1 - sin(((i - 0.5) / napodize) * pi * 0.5);
%     % top
%     image(i, :) = image(i, :) + diff1 .* fact;
%     % bottom
%     image(ny + 1 - i, :) = image(ny + 1 - i, :) + diff2 .* fact;
% end
% 
% % Column-wise soften
% diff1 = 0 - image(:, 1);      % Difference between 0 and the first (1) column;
% diff2 = 0 - image(:, nx);     % Difference between 0 and the last (512) column;
% for j = 1: napodize
%     fact = 1 - sin(((j - 0.5) / napodize) * pi * 0.5);
%     % left
%     image(:, j) = image(:, j) + diff1 .* fact;
%     % right
%     image(:, nx + 1 - j) = image(:, nx + 1 - j) + diff2 .* fact;
% end
% 
% image(image < 0) = 0;



% My algorithm # 2

fact = 1 / napodize * pi * 0.5;
% Row-wise soften
for i = 1: napodize
    % top
    image(i, :) = image(i, :) * (sin(fact * (i - 1))) .^ 2;
    % bottom
    image(ny + 1 - i, :) = image(ny + 1 - i, :) * (sin(fact * (i - 1))) .^ 2;
end

% Column-wise soften
for j = 1: napodize
    % left
    image(:, j) = image(:, j) * (sin(fact * (j - 1))) .^ 2;
    % right
    image(:, nx + 1 - j) = image(:, nx + 1 - j) * (sin(fact * (j - 1))) .^ 2;
end

image(image < 0) = 0;

end

%% to apodize (or edge softening) every 2D slice
function image = cosapodize(nx, ny, image)

x = 1 : nx;
y = 1 : ny;
[xfact, yfact] = meshgrid(x, y);
xfact = single(sin(pi * (xfact - 0.5) / nx));
yfact = single(sin(pi * (yfact - 0.5) / ny));
image = image .* xfact .* yfact;

image(image < 0) = 0;

end

%% Called by SIM_Reconstructor::loadAndRescaleImage() to do bleach correction for obj.m_reconData.savedBands(: ,:, :, phase, direction)
function data = rescaleDriver(it, iw, zoffset, params, imgParams, data)
%  data assumed taken with z changing fast, direction changing slowly

for direction = 1 : params.ndirs
    rawImages = data.savedBands(:, :, :, :, direction);
    
    for z = 1 : imgParams.nz
        sum_buffer = zeros([params.nphases, 1], 'single');
        for phase = 1 : params.nphases
            % Get sum value of current z image for individual 5 phases
            sum_buffer(phase) = sum(rawImages(:, :, z + zoffset, phase), 'all');
        end
        
        if (direction == 1) && ~(params.equalizet && it ~= 0)
            % That is to say, sum_dir0_phase0 will save the sum value of
            % individual z image for direction = 1 and phase = 1 when t = 0
            data.sum_dir0_phase0(iw * imgParams.nz + z) = sum_buffer(1);
        end
        
%         if (params.equalizez) && ~(params.bTwolens)
        if (params.equalizez)
            % sum value of image intensity under direction = 1, phase = 1, z = 1
            ref = data.sum_dir0_phase0(iw * imgParams.nz + 1);
        else
            % sum value of image intensity under direction = 1, phase = 1
            ref = data.sum_dir0_phase0(iw * imgParams.nz + z);
        end
        
        for phase = 1 : params.nphases
            ratio = ref / sum_buffer(phase);
            data.savedBands(:, :, z + zoffset, phase, direction) = data.savedBands(:, :, z + zoffset, phase, direction) .* ratio;
        end
    end
end

end

%% Save wide field image
function saveWidefield(it, iw, params, imgParams, images)

wide_field = sum(images, [4, 5]) / params.nphases / params.ndirs;
[Height, Width, Depth] = size(wide_field);
if (Depth == 1)
    wide_field = imresize(wide_field, [Height * params.zoomfact, Width * params.zoomfact]);
else
    wide_field = imresize3(wide_field, [Height * params.zoomfact, Width * params.zoomfact, Depth * params.z_zoom]);
end

wide_field_fft = log(abs(fftshift(fftn(wide_field))));

% Write wide-field image
fileID = fopen(strcat(params.ifiles, '\wide_field_', num2str(it), '.tif'), 'w+');
[Height, Width, Depth] = size(wide_field);
minval = min(wide_field, [], 'all');
maxval = max(wide_field, [], 'all');
header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(wide_field), 1, Depth, 1, imgParams.dz / params.z_zoom, minval, maxval, imgParams.dy / params.zoomfact);
% header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(wide_field), 1, Depth, 1, imgParams.dz, minval, maxval, imgParams.dy);
wide_field = permute(wide_field, [2 1 3]);
wide_field = reshape(wide_field, [1, Height * Width * Depth]);

[~, ~, system_endian] = computer;
if system_endian ~= header.endian
    wide_field = swapbytes(wide_field);
end

wide_field = typecast(wide_field, 'uint8');
fwrite(fileID, wide_field, 'uint8');
fclose(fileID);

% Wrtie wide-field fft data image
if (params.bsavefft)
    fileID = fopen(strcat('wide_field_fft_', num2str(it), '.tif'), 'w+');
    [Height, Width, Depth] = size(wide_field_fft);
    minval = min(wide_field_fft, [], 'all');
    maxval = max(wide_field_fft, [], 'all');
%     header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(wide_field_fft), 1, Depth, 1, imgParams.dz / params.z_zoom, minval, maxval, imgParams.dy / params.zoomfact);
    header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(wide_field_fft), 1, Depth, 1, imgParams.dz, minval, maxval, imgParams.dy);
    wide_field_fft = permute(wide_field_fft, [2 1 3]);
    wide_field_fft = reshape(wide_field_fft, [1, Height * Width * Depth]);
    
    if system_endian ~= header.endian
        wide_field_fft = swapbytes(wide_field_fft);
    end
    
    wide_field_fft = typecast(wide_field_fft, 'uint8');
    fwrite(fileID, wide_field_fft, 'uint8');
    fclose(fileID);
end

end

%% Save Richardson-Lucy deconvolved wide field image
function saveRLWidefield(it, iw, otf, params, imgParams, images)

wide_field = sum(images, [4, 5]) / params.nphases / params.ndirs;
outputOtf = complex(zeros(size(wide_field), 'single'), zeros(size(wide_field), 'single'));

% Calculate krscale
dkr = 1.0 / (imgParams.ny * imgParams.dy);
krscale = dkr / otf.dkrotf;

% Calculate kzscale
if (imgParams.dz > 0)
    dkz = 1.0 / (imgParams.nz0 * imgParams.dz);
else
    dkz = otf.dkzotf;
end
% ratio of axial direction pixel scales of data and otf
kzscale = dkz / otf.dkzotf;

lambdaem = (imgParams.wave / params.nimm) / 1.0e3;
alpha = asin(params.na / params.nimm);
zdistcutoff = ceil(((1 - cos(alpha)) / lambdaem) / dkz);

if (zdistcutoff >= fix(imgParams.nz0 / 2))
    if ((fix(imgParams.nz0 / 2) - 1) > 0)
        zdistcutoff = fix(imgParams.nz0 / 2) - 1;   % 3D-SIM
    else
        zdistcutoff = 0;                            % 2D-SIM
    end
end

otfval = otf.getOtfMatrix3D(0, 1, imgParams.nx, imgParams.ny, 0, 0, krscale, zdistcutoff, kzscale);

z0 = -zdistcutoff : zdistcutoff;
iz = mod((z0 + imgParams.nz0), imgParams.nz0); 
outputOtf(:, :, iz + 1) = otfval;


if (true)
    outputOtf_buffer = log(abs(fftshift(outputOtf)));
    fileID = fopen('wide_field_RL_OTF.tif', 'w+');
    [Height, Width, Depth] = size(outputOtf_buffer);
    minval = min(outputOtf_buffer, [], 'all');
    maxval = max(outputOtf_buffer, [], 'all');
    header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(outputOtf_buffer), 1, Depth, 1, imgParams.dz, minval, maxval, imgParams.dy);
    outputOtf_buffer = permute(outputOtf_buffer, [2 1 3]);
    outputOtf_buffer = reshape(outputOtf_buffer, [1, Height * Width * Depth]);
    
    [~, ~, system_endian] = computer;
    if system_endian ~= header.endian
        outputOtf_buffer = swapbytes(outputOtf_buffer);
    end
    
    outputOtf_buffer = typecast(outputOtf_buffer, 'uint8');
    fwrite(fileID, outputOtf_buffer, 'uint8');
    fclose(fileID);
end


% Run Richardson-Lucy deconvolution
wide_field_RL = deconvolve(wide_field, outputOtf, params.iter_input, false);
wide_field_RL_fft = log(abs(fftshift(fftn(wide_field_RL))));

% Write RL-deconvolved wide-field image
fileID = fopen(strcat('wide_field_RL_', num2str(it), '.tif'), 'w+');
[Height, Width, Depth] = size(wide_field_RL);
minval = min(wide_field_RL, [], 'all');
maxval = max(wide_field_RL, [], 'all');
header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(wide_field_RL), 1, Depth, 1, imgParams.dz, minval, maxval, imgParams.dy);
wide_field_RL = permute(wide_field_RL, [2 1 3]);
wide_field_RL = reshape(wide_field_RL, [1, Height * Width * Depth]);

[~, ~, system_endian] = computer;
if system_endian ~= header.endian
    wide_field_RL = swapbytes(wide_field_RL);
end

wide_field_RL = typecast(wide_field_RL, 'uint8');
fwrite(fileID, wide_field_RL, 'uint8');
fclose(fileID);

% Wrtie RL-deconvolved wide-field fft data image
if (params.bsavefft)
    fileID = fopen(strcat('wide_field_RL_fft_', num2str(it), '.tif'), 'w+');
    [Height, Width, Depth] = size(wide_field_RL_fft);
    minval = min(wide_field_RL_fft, [], 'all');
    maxval = max(wide_field_RL_fft, [], 'all');
    header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(wide_field_RL_fft), 1, Depth, 1, imgParams.dz, minval, maxval, imgParams.dy);
    wide_field_RL_fft = permute(wide_field_RL_fft, [2 1 3]);
    wide_field_RL_fft = reshape(wide_field_RL_fft, [1, Height * Width * Depth]);
    
    if system_endian ~= header.endian
        wide_field_RL_fft = swapbytes(wide_field_RL_fft);
    end
    
    wide_field_RL_fft = typecast(wide_field_RL_fft, 'uint8');
    fwrite(fileID, wide_field_RL_fft, 'uint8');
    fclose(fileID);
end

end

%% Run 2D Richardson-Lucy deconvolution on each input image
function bands = RL_inputdriver(otf, zoffset, params, imgParams, bands)

dkr = 1.0 / (imgParams.ny * imgParams.dy);
krscale = dkr / otf.dkrotf;
inputOtf = otf.getOtfMatrix2D(0, 1, imgParams.nx, imgParams.ny, 0, 0, krscale);

if (true)
    inputOtf_buffer = log(abs(fftshift(inputOtf)));
    fileID = fopen('input_RL_2D_OTF.tif', 'w+');
    [Height, Width] = size(inputOtf_buffer);
    minval = min(inputOtf_buffer, [], 'all');
    maxval = max(inputOtf_buffer, [], 'all');
    header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(inputOtf_buffer), 1, 1, 1, imgParams.dz, minval, maxval, imgParams.dy);
    inputOtf_buffer = permute(inputOtf_buffer, [2 1]);
    inputOtf_buffer = reshape(inputOtf_buffer, [1, Height * Width]);
    
    [~, ~, system_endian] = computer;
    if system_endian ~= header.endian
        inputOtf_buffer = swapbytes(inputOtf_buffer);
    end
    
    inputOtf_buffer = typecast(inputOtf_buffer, 'uint8');
    fwrite(fileID, inputOtf_buffer, 'uint8');
    fclose(fileID);
end

for direction = 1 : params.ndirs
    for z = 1 : imgParams.nz
        for phase = 1 : params.nphases
                bands(:, :, z + zoffset, phase, direction) = deconvolve(bands(:, :, z + zoffset, phase, direction), inputOtf, params.iter_input, false);
        end
    end
end

end

%% Save Richardson-Lucy deconvolution input images
function saveRLinput(it, iw, params, imgParams, images)

minval = min(images, [], 'all');
maxval = max(images, [], 'all');
images = permute(images, [2 1 4 5 3]);
images = images(:);
images = reshape(images, [1, size(images)]);

fileID = fopen(strcat('RL_input_', num2str(it), '.tif'), 'w+');
header = ImageJ_formatted_TIFF.write_IFD(fileID, imgParams.nx, imgParams.ny, class(images), ...
    params.nphases * params.ndirs , imgParams.nz0, 1, imgParams.dz, minval, maxval, imgParams.dy);

[~, ~, system_endian] = computer;
if system_endian ~= header.endian
    images = swapbytes(images);
end

images = typecast(images, 'uint8');
fwrite(fileID, images, 'uint8');
fclose(fileID);

end





%% Find the modulation vectors, modulation amplitudes and phases for all directions
function [params, driftParams, data] = findModulationVectorsAndPhasesForAllDirections(zoffset, otf, params, imgParams, driftParams, data)

% 2D in-place FFT of every 2D slice
data = transformXYSlice(zoffset, params, imgParams, data);
% After this function, m_reconData.savedBands[direction][phase] has been 2D FFT for nz layers;
% In CUDA: m_reconData.savedBands[direction][phase]((514 * 512 * float)* nz0) --> m_reconData.savedBands[direction][phase]((512 * 257 * FloatComplex)* nz0);

% loop through pattern directions
for direction = 1 : params.ndirs
    rawImages = data.savedBands(:, :, :, :, direction);
    
    phaseList = zeros([params.nphases, 1], 'single');
    if (params.phaseSteps ~= 0)
        % User specified the non-ideal (i.e., not 2pi/5) phase steps for each orientation.
        % Under this circumstance, calculate the non-ideal sepMatrix:
        for i = 1 : params.nphases
            phaseList(i) = (i - 1) * params.phaseSteps(direction);
        end
        [data.sepMatrix, data.noiseVarFactors] = makematrix(params.nphases, params.norders, direction, phaseList, data.sepMatrix, data.noiseVarFactors);
    end
    
    % Unmixing info components in real or reciprocal space:
    data.savedBands(:, :, :, :, direction) = separate(params.nphases, params.norders, rawImages, data.sepMatrix);
    % After this function, m_reconData.savedBands[direction][order = 0, 1re, 1im, 2re, 2im] save the unmixed lateral (2D) frequency components,
    % which are defined by the matrix multiplication result between "m_reconData.sepMatrix" and "m_reconData.savedBands[direction]"
    
    fprintf("Direction %d:\n", direction);
    
    if (imgParams.nz > 1)
        % 1D FFT of a stack of 2D FFTs along the 3rd dimension to obtain equivalent of 3D FFT
        for i = 1 : params.nphases
            data.savedBands(:, :, :, i, direction) = fft(data.savedBands(:, :, :, i, direction), imgParams.nz0, 3);
        end
        % After these 5 phases loops, m_reconData.savedBands[direction][order = 0, 1re, 1im, 2re, 2im] save the unmixed 3D frequency components;
    end
    
    bands = data.savedBands(:, :, :, :, direction);
    
    % Use the OTF to simulate an ideal point source; replace bands with simulated data
    if (params.bMakemodel)
        makemodeldata(data.savedBands(:, :, :, :, direction), data.k0(direction), otf, params, imgParams);
    end
    
    % save the separated bands (in frequency domain) if requested
    if (params.bSaveSeparated)
        fprintf("Saving separated bands\n");
        saveSeparatedbands(params, imgParams, direction, bands);
        fprintf("Separated bands saved\n\n");
        %         continue;   % skip the k0 search and modamp fitting
    end
    
    %{
    After separation and FFT, rawImages, now referred to as bands,
    * contains the center band (bands[0]),
    * the real (bands[1], bands[3], etc.) and
    * imaginary part (bands[2], bands[4], etc.) bands
    * of the different orders
    %}
    
    
    % If we do not know anything about k0guess
    if (~(imgParams.ntimes > 1 && imgParams.curTimeIdx > 0 && params.bUseTime0k0))
        if (isempty(data.k0guess(direction).x))
            fprintf("Calculating k0guess(direction %d):\n", direction);
            data.k0guess(direction) = estimatek0guess(bands, otf, params, imgParams);
        end
        data.k0(direction) = data.k0guess(direction);
        fprintf("k0guess(direction %d) = (%f, %f)\n\n", direction, data.k0guess(direction).x, data.k0guess(direction).y);
    end
    
    % Now to fix 3D drift between dirs estimated by determinedrift_3D()
    if (direction ~= 1 && params.bFixdrift)
        fixdrift_bt_dirs(bands, driftParams.drift_bt_dirs(direction), imgParams);
    end
    
    % assume k0 vector not well known, so fit for it
    amp_inv = complex(single(0.0), single(0.0));
    amp_combo = complex(single(0.0), single(0.0));
    
    % dir_ is used in upcoming calls involving data->otf, to differentiate the cases of common OTF and dir-specific OTF
    dir_ = 1;   % Suppose common OTF for all angles / directions
    if (params.bOneOTFperAngle)
        dir_ = direction;
    end
    
    if (params.bSearchforvector && ~(imgParams.ntimes > 1 && imgParams.curTimeIdx > 0 && params.bUseTime0k0))
        % In time series, can choose to use the time 0 k0 fit for the rest of the series.
        % Find initial estimate of modulation wave vector m_reconData.k0(direction) by cross-correlation.
        [data.overlap0, data.overlap1, data.k0(direction)] = findk0(bands, data.overlap0, data.overlap1, ...
            data.k0(direction), otf, dir_, params, imgParams);
        
        fprintf("Initial guess by findk0() of k0[direction %d] = (%f, %f)\n\n", direction, data.k0(direction).x, data.k0(direction).y);
        
        % refine the cross-corr estimate of k0 by a search for best k0
        % vector direction and magnitude using real space waves
        
        fprintf("before fitk0andmodamp\n");
        params.recalcarrays = 3;
        % if k0 is very close to the guess, we can save time by not recalculating the overlap arrays
        
        % (1): to modify / update m_reconData.k0[direction].x and y
        % (2): to compute m_reconData.amp[direction][order = 1] and [order = 2];
        [data.k0(direction), data.amp(:, direction), data.overlap0, data.overlap1] = fitk0andmodamps(bands, data.overlap0, data.overlap1, ...
            data.k0(direction), data.amp(:, direction), otf, dir_, direction, params, imgParams);
        
        if imgParams.curTimeIdx == 0
            data.k0_time0(direction) = data.k0(direction);
        end
        
        % check if the k0 vector found is reasonably close to the guess
        K0_WARNING_THRESH = 2;
        deltak0 = struct('x',[], 'y', []);
        deltak0.x = data.k0(direction).x - data.k0guess(direction).x;
        deltak0.y = data.k0(direction).y - data.k0guess(direction).y;
        dist = sqrt(deltak0.x ^ 2 + deltak0.y ^ 2);
        if (dist > K0_WARNING_THRESH)
            warning("best fit for k0 is %f pixels from expected value.\n", dist);
        end
        
        if (imgParams.ntimes > 1 && imgParams.curTimeIdx > 0 && dist > 2 * K0_WARNING_THRESH)
            data.k0(direction) = data.k0_time0(direction);
            fprintf("k0 estimate of time point 0 is used instead\n");
            
            for order = 1 : (params.norders - 1)
                if (imgParams.nz0 > 1)
                    [data.amp(order + 1, direction), ~, ~, corr_coeff, data.overlap0, data.overlap1] = findrealspacemodamp(bands, data.overlap0, data.overlap1, ...
                        0, order, data.k0(direction), true, otf, dir_, params, imgParams);
                else
                    [data.amp(order + 1, direction), ~, ~, corr_coeff, data.overlap0, data.overlap1] = findrealspacemodamp(bands, data.overlap0, data.overlap1, ...
                        order - 1, order, data.k0(direction), true, otf, dir_, params, imgParams);
                end
                
                fprintf("modamp mag=%f, phase=%f\n, correlation coeff=%f\n\n", abs(data.amp(order + 1, direction)), angle(data.amp(order + 1, direction)), corr_coeff);                
            end
        end
        
    else
        % assume k0 vector well known, so just fit for the modulation amplitude and phase
        fprintf("known k0 for direction %d = (%f, %f) \n", direction, data.k0(direction).x, data.k0(direction).y);
        
        for order = 1 : (params.norders - 1)
            [data.amp(order + 1, direction), amp_inv, amp_combo, corr_coeff, data.overlap0, data.overlap1] = findrealspacemodamp(bands, data.overlap0, data.overlap1, ...
                0, order, data.k0(direction), true, otf, dir_, params, imgParams);
            fprintf("modamp mag=%f, phase=%f\n", abs(data.amp(order + 1, direction)), angle(data.amp(order + 1, direction)));
            fprintf("reverse modamp mag=%f, phase=%f\n", 1.0 / abs(amp_inv), -angle(amp_inv));
            fprintf("combined modamp mag=%f, phase=%f\n", abs(amp_combo), angle(amp_combo));
            fprintf("correlation coeff=%f\n\n", corr_coeff);
            
            if (params.bSaveOverlaps)
                fprintf("Saving overlapped bands\n");
                saveOverlappedbands(params, imgParams, direction, order, data.overlap0, data.overlap1);
                fprintf("Overlapped bands saved\n\n");
            end            
        end
    end
    
    % Not so useful in 3D SIM case.
    if (imgParams.nz == 1)
        % In 2D SIM, amp stores modamp's between each adjacent pair of bands.
        % We want to convert this to modamp w.r.t. order 0
        for order = 2 : (params.norders - 1)
            data.amp(order + 1, direction) = data.amp(order + 1, direction) * data.amp(order, direction);
        end
    end
    
    % Also not so useful in real application, we prefer to use modamp's amplitude based on calculation.
    if ~isempty(params.forceamp)
        % force modamp's amplitude to be a value user provided (ideally should be 1)
        for order = 1 : (params.norders - 1)
            a = abs(data.amp(order + 1, direction));
            if a < params.forceamp(order)
                ampfact = params.forceamp(order) / a;
                data.amp(order + 1, direction) = data.amp(order + 1, direction) * ampfact;
                fprintf("modamp mag=%f, phase=%f  \n", abs(data.amp(order + 1, direction)), angle(data.amp(order + 1, direction)));
            end
        end
    end
    
    % In 2D NLSIM, we often don't trust the modamp fit between neighboring high-order components;
    % only if fitallphases is True do all fitted modamps get actually used;
    % otherwise, order 2 and above's global phase is inferred from order 1's phase.
    if (~params.bFitallphases)
        base_phase = angle(data.amp(2, direction));     % base_phase is calculated from m_reconData.amp[direction][order = 1];
        for order = 2 : (params.norders - 1)
            phi = order * base_phase;
            data.amp(order + 1, direction) = abs(data.amp(order + 1, direction)) * exp(1j * phi);
        end
    end
end

end

%% to apply 2D FFT of every 2D slice
function data = transformXYSlice(zoffset, params, imgParams, data)
% In CUDA: m_reconData.savedBands[direction][phase]((514 * 512 * float) * nz0) --> m_reconData.savedBands[direction][phase]((257 * 512 * cuFloatComplex) * nz0);

for direction = 1 : params.ndirs
    for phase = 1 : params.nphases
        for z = 1 : imgParams.nz
            data.savedBands(:, :, z + zoffset, phase, direction) = fft2(data.savedBands(:, :, z + zoffset, phase, direction));
        end
    end
end

end

%% to unmix info components in real or reciprocal space
function output = separate(nphases, norders, rawImages, sepMatrix)
%{
Noted that rawImages = m_reconData.savedBands[direction] here has already been the frequency domain; (2D-FFT for each z layer already)

				phase = 0;		2pi / 5					2 * 2pi / 5				3 * 2pi / 5				4 * 2pi / 5
m_reconData.sepMatrix = [1,			1,						1,						1,						1;						order = 0
                        1, cos(1 * 1 * 2pi / 5),   cos(2 * 1 * 2pi / 5),   cos(3 * 1 * 2pi / 5),	cos(4 * 1 * 2pi / 5);			order = 1, bandre
						0, sin(1 * 1 * 2pi / 5),   sin(2 * 1 * 2pi / 5),   sin(3 * 1 * 2pi / 5),	sin(4 * 2 * 2pi / 5);			order = 1, bandim
						1, cos(1 * 2 * 2pi / 5),   cos(2 * 2 * 2pi / 5),   cos(3 * 2 * 2pi / 5),	cos(4 * 2 * 2pi / 5);			order = 2, bandre
						0, sin(1 * 2 * 2pi / 5),   sin(2 * 2 * 2pi / 5),   sin(3 * 2 * 2pi / 5),	sin(4 * 2 * 2pi / 5);]			order = 2, bandim

%}
output = complex(zeros(size(rawImages), 'single'), zeros(size(rawImages), 'single'));

for i = 1 : (norders * 2 - 1)
    for j = 1 : nphases
        output(:, :, :, i) = output(:, :, :, i) + rawImages(:, :, :, j) .* sepMatrix(i, j);
    end
end

%{
    After the 5 loops for orders:
	output[order = 0][offset] = m_reconData.sepMatrix[0][1][2][3][4] * m_reconData.savedBands[direction][phase: 0,1,2,3,4][offset];
	output[order = 1re][offset] = m_reconData.sepMatrix[5][6][7][8][9] * m_reconData.savedBands[direction][phase: 0,1,2,3,4][offset];
	output[order = 1im][offset] = m_reconData.sepMatrix[10][11][12][13][14] * m_reconData.savedBands[direction][phase: 0,1,2,3,4][offset];
	output[order = 2re][offset] = m_reconData.sepMatrix[15][16][17][18][19] * m_reconData.savedBands[direction][phase: 0,1,2,3,4][offset];
	output[order = 2im][offset] = m_reconData.sepMatrix[20][21][22][23][24] * m_reconData.savedBands[direction][phase: 0,1,2,3,4][offset];
%}

end

%% to use the OTF to simulate an ideal point source; replace bands with simulated data
function makemodeldata(bands, k0, otf, params, imgParams)
fprintf("In makemodeldata.\n");
end

%% to fix 3D drift between dirs estimated by determinedrift_3D()
function fixdrift_bt_dirs(bands, drift, imgParams)
	fprintf("In fixdrift_bt_dirs.\n");
end

%% Save separated bands fft data image
function saveSeparatedbands(params, imgParams, direction, images)

separated_bands = zeros([imgParams.ny, imgParams.nx, imgParams.nz0, 2 * params.norders - 1], 'single');
for order = 0 : (params.norders - 1)
    if order == 0
        separated_bands(:, :, :, order + 1) = log(abs(fftshift(images(:, :, :, order + 1))));
    else
        bandre = fftshift(images(:, :, :, 2 * order));
        bandim = fftshift(images(:, :, :, 2 * order + 1));
        separated_bands(:, :, :, 2 * order) = log(abs(bandre + 1j * bandim));       % Minus bands
        separated_bands(:, :, :, 2 * order + 1) = log(abs(bandre - 1j * bandim));   % Plus bands
    end
end

fileID = fopen(strcat(params.fileSeparated, '_dir_', num2str(direction), '.tif'), 'w+');
[Height, Width, Depth, Channel] = size(separated_bands);
minval = min(separated_bands, [], 'all');
maxval = max(separated_bands, [], 'all');
header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(separated_bands), Channel, Depth, 1, imgParams.dz, minval, maxval, imgParams.dy);
separated_bands = permute(separated_bands, [2 1 4 3]);
separated_bands = reshape(separated_bands, [1, Height * Width * Channel * Depth]);

[~, ~, system_endian] = computer;
if system_endian ~= header.endian
    separated_bands = swapbytes(separated_bands);
end

separated_bands = typecast(separated_bands, 'uint8');
fwrite(fileID, separated_bands, 'uint8');
fclose(fileID);

end

%% to estimate k0guess if we know nothing about it
function k0guess = estimatek0guess(bands, otf, params, imgParams)
%{
    bands = m_reconData.savedBands(:, :, :, :, direction)
    otf = obj.otf;
    params = m_myParams;
%}

fitorder1 = 0;
if (imgParams.nz > 1)
    fitorder2 = 2;  % 3D-SIM
else
    fitorder2 = 1;  % 2D-SIM
end

% Calculate krscale
dkr = 1.0 / (imgParams.ny * imgParams.dy);
krscale = dkr / otf.dkrotf;
rdistcutoff = fix((params.na * 2.0 / (imgParams.wave / 1.0e3)) / dkr);

% Calculate kzscale
if (imgParams.dz > 0)
    dkz = 1.0 / (imgParams.nz0 * imgParams.dz);
else
    dkz = otf.dkzotf;
end

lambdaem = (imgParams.wave / params.nimm) / 1.0e3;
alpha = asin(params.na / params.nimm);
zdistcutoff = ceil(((1 - cos(alpha)) / lambdaem) / dkz);

if (zdistcutoff >= fix(imgParams.nz0 / 2))
    if ((fix(imgParams.nz0 / 2) - 1) > 0)
        zdistcutoff = fix(imgParams.nz0 / 2) - 1;   % 3D-SIM
    else
        zdistcutoff = 0;                            % 2D-SIM
    end
end
z0 = -zdistcutoff : zdistcutoff;
iz = mod((z0 + imgParams.nz0), imgParams.nz0);

fprintf("fitorder2=%d, rdistcutoff=%d, zdistcutoff=%f\n", fitorder2, rdistcutoff, zdistcutoff);

% c1: band(order = 0); c2: band(order = +1 or +2)
c1 = bands(:, :, iz + 1, fitorder1 + 1);
c2 = bands(:, :, iz + 1, fitorder2 * 2) - 1j .* bands(:, :, iz + 1, fitorder2 * 2 + 1);
c1 = ifftn(c1);
c2 = ifftn(c2);
crosscorr_c = sum(c2 .* conj(c1), 3);
crosscorr_c = fft2(crosscorr_c);

if (imgParams.nz > 1)
    minDist = 0.6 * otf.getCutoff() / dkr;  % 3D-SIM
else
    minDist = 0.2 * otf.getCutoff() / dkr;  % 2D-SIM
end

k0guess = struct('x',[], 'y', []);
[k0guess.x, k0guess.y, ~, ~] = locatePeak(crosscorr_c, minDist);

k0guess.x = k0guess.x / fitorder2;
k0guess.y = k0guess.y / fitorder2;

end

%% Locates position, magnitute and phase of the highest peak in 'mat'
function [xPos, yPos, big, phase] = locatePeak(mat, kMin)

ny = size(mat, 1);
nx = size(mat, 2);

x = 0 : (nx - 1);
y = 0 : (ny - 1);
x(x >= (nx / 2)) = x(x >= (nx / 2)) - nx;
y(y >= (ny / 2)) = y(y >= (ny / 2)) - ny;
[X, Y] = meshgrid(x, y);
rad = sqrt(X .* X + Y .* Y);

mat(rad <= kMin) = 0;
big = max(abs(mat), [], 'all');
[yPos, xPos] = find(abs(mat) == big);
phase = angle(mat(yPos, xPos));

yPos = yPos - 1;
xPos = xPos - 1;

if (xPos > fix(nx / 2))
    xPos = xPos - nx;
end

if (yPos > fix(ny / 2))
    yPos = yPos - ny;
end

end

%% to find the initial estimate of modulation wave vector m_reconData.k0[direction] by cross-correlation */
function [overlap0, overlap1, k0] = findk0(bands, overlap0, overlap1, k0, otf, direction, params, imgParams)
%{
    bands = m_reconData.savedBands(:, :, :, :, direction)
    overlap0 = class: 'single complex'; size:[m_imgParams.ny, m_imgParams.nx, m_imgParams.nz]
    otf = obj.otf;
    direction = dir_ = 1;
    params = obj.m_myParams;
    imgParams = obj.m_imgParams;
    k0 = m_reconData.k0(direction); Suppose direction = 1, k0.x = m_reconData.k0guess(1).x = 100.438431; k0.y = -0.803525
%}

fitorder1 = 0;
if (imgParams.nz0 > 1)
    fitorder2 = 2;  % 3D-SIM
else
    fitorder2 = 1;  % 2D-SIM
end

% (1) calculate overlap0 = "m_reconData.savedBands[direction][order = 0] multiplied with order2 otf shifted by (kx, ky)"; that is: D0(k)OTF2(k + 2P);
% (2) calculate overlap1 = "m_reconData.savedBands[direction][order = +2] multiplied with order0 otf shifted by (-kx, -ky)"; that is: D+2(k)OTF0(k - 2P);
% (3) apply in-place ifft transform of overlap0 and overlap1;
 [overlap0, overlap1] = makeoverlaps(bands, overlap0, overlap1, fitorder1, fitorder2, k0.x, k0.y, otf, direction, params, imgParams);
%{
Cross-correlation Defination & Properties: 
	1. The cross-correlation of functions f(t) and g(t) is equivalent to the convolution of conjugate of f(-t) and g(t);
	2. Analogous to the convolution theorem, the cross-correlation satisfies: 
				F[f(t) cross* g(t)] = F[conj[f(-t)] conv* g(t)] = F[conj[f(-t)]] * F[g(t)]] = conj(F[f(t)]) * F[g(t)];
				
	In our application, we want to calculate the cross-correlation of overlap0 & overlap1, which are all in Fourier domain; 
	So all the Fourier transform in 2 should be modified by inverse Fourier transform:
				ifft[overlap1 cross* overlap0] = conj(ifft[overlap1]) * ifft[overlap0];
	That is the reason why we apply in-place ifft transform of overlap0 and overlap1 here.
	Therefore, [overlap1 cross* overlap0] = fft[conj(ifft[overlap1]) * ifft[overlap0]];
	
	After in-place ifft transform of overlap0 and overlap1, cross_correlation = fft[overlap1 * conj(overlap0)];
	
	Noted that: The cross-correlation function is not symmetric, xcorr(x1, x2) != xcorr(x2, x1) !!!
				Therefore, we should pay attention to the frequency domain calculation to see who takes the conjuga.
				Simply remember that the signal is used as a reference, which one is conjugated.
	In our case, we want overlap0 as the reference, and "shift" overlap1.
%}

% crosscorr_c is the complex form of cross-correlation between overlap0 and overlap1;
% Decrease the dimension from 3D to 2D, add all the data in z-direction;
crosscorr_c = sum(overlap1 .* conj(overlap0), 3);	

% In-place fft transform of crosscorr_c;
crosscorr_c = fft2(crosscorr_c);

old_k0 = k0;
intensities = abs(crosscorr_c) .^ 2;
k0 = findpeak(intensities, imgParams.nx, imgParams.ny, k0);
% k0.x, k0.y are within range (-1, 512)

if ((k0.x - old_k0.x) > fix(imgParams.nx / 2))
    k0.x = k0.x - imgParams.nx;
end
if ((old_k0.x - k0.x) > fix(imgParams.nx / 2))
    k0.x = k0.x + imgParams.nx;
end

if ((k0.y - old_k0.y) > fix(imgParams.ny / 2))
    k0.y = k0.y - imgParams.ny;
end
if ((old_k0.y - k0.y) > fix(imgParams.ny / 2))
    k0.y = k0.y + imgParams.ny;
end

% return k0 of the first order, no matter which fitorder2 is used */
k0.x = k0.x / fitorder2;
k0.y = k0.y / fitorder2;

end

%% to calculate overlap0 = ifft[D0(k)OTF2(k + 2P)], overlap1 = ifft[D+2(k)OTF0(k - 2P)] & overlap1 = ifft[D+2(k)OTF0(k - P)]
function [overlap0, overlap1] = makeoverlaps(bands, overlap0, overlap1, order1, order2, k0x, k0y, otf, direction, params, imgParams)
% (1) calculate "m_reconData.savedBands[direction][order = 0] multiplied with order2 otf shifted by (kx, ky)", and saved into overlap0;
% (2) calculate "m_reconData.savedBands[direction][order = +2] multiplied with order0 otf shifted by (-kx, -ky)", and saved into overlap1; 
% (3) apply in-place ifft transform of overlap0 and overlap1;

%{
    bands = obj.m_reconData.savedBands(:, :, :, :, direction)
    overlap0 = class: 'single complex'; size:[m_imgParams.ny, m_imgParams.nx, m_imgParams.nz]
    overlap1 = class: 'single complex'; size:[m_imgParams.ny, m_imgParams.nx, m_imgParams.nz]
    order1 = fitorder1 = 0;
	order2 = fitorder2 = 2; 
    (or order2 = 1 to compute order = 1 modulation factor;)
    k0x = m_reconData.k0guess(1).x = 100.438431;
    k0y = m_reconData.k0guess(1).y = -0.803525;
    otf = obj.otf;
    direction = dir_ = 1;
    params = obj.m_myParams;
    imgParams = obj.m_imgParams;
%}


order0_2_factor = 1.0;
if (imgParams.nz0 > 1)
    order0_2_factor = 5.0;
end


% "rdistcutoff" is the lateral pixels of maximum frequency of emmision wavelength in current matrix.
rdistcutoff = fix((params.na * 2.0 / (imgParams.wave(1) / 1.0e3)) / imgParams.dkr);
%{
rdistcutoff = (m_myParams.na * 2.0 / (525 nm / 1.0e3)) / (1 / (40.96 um));
			= (1.35 * 2.0 / 0.525 um) * 40.96 um = 210 pixels;
"2NA / wave" is the abbe limit of the objective in frequency space for emmision; "1 / 2dy" is the maximum lateral frequency in the matrix;
rdistcutoff = 2NA / wave / (1 / 2dy) * ny_half = 2NA / wave / (1 / (dy * ny)) = 2NA / wave / dkr;
%}

if (rdistcutoff > fix(imgParams.nx / 2))
    rdistcutoff = fix(imgParams.nx / 2);                % rdistcutoff = 210 < 256 = nx / 2;
end

k0pix = sqrt(k0x ^ 2 + k0y ^ 2);
k0pix = k0pix * (params.norders - 1);
% k0 magnitude (for highest order) in inverse microns
k0mag = k0pix * imgParams.dkr;
lambdaem = (imgParams.wave(1) / params.nimm) / 1.0e3;	% lambdaem = (wave / m_myParams.nimm) / 1.0e3 = (525 nm / 1.406) / 1.0e3 = (0.525 um / 1.406);
lambdaexc = (params.lambda / params.nimm) / 1.0e3;      % lambdaexc = (488 nm / 1.406) / 1.0e3 = (0.488 um / 1.406);
% lambdaexc_medium is to determine standing wave position in frequency domain
lambdaexc_medium = params.lambda / params.nmed / 1.0e3; % lambdaexc_medium = (488 nm / 1.333) / 1.0e3 = (0.488 um / 1.333);
% alpha is the aperture angle of objectives
alpha = asin(params.na / params.nimm);                  % alpha = asin(1.35 / 1.406) = 1.2876 radians = 73.7748 degrees;
% beta is the angle of center of side illumination beams
beta = asin(k0mag / (2.0 / lambdaexc_medium));



% "zdistcutoff" is the axial pixels of maximum frequency of emmision OTF;
if ~(params.bTwolens)
    zdistcutoff = ceil(((1.0 - cos(alpha)) / lambdaem) / imgParams.dkz);
else
    if (order2 == 1)
        zdistcutoff = ceil(((1.0 - cos(alpha)) / lambdaem + (1 + cos(beta)) / lambdaexc) / imgParams.dkz);
%         zdistcutoff = ceil(((1.0 - cos(alpha)) / lambdaem) / imgParams.dkz);
    elseif (order2 == 2)
        zdistcutoff = ceil(((1.0 - cos(alpha)) / lambdaem) / imgParams.dkz);
    end
end

if (zdistcutoff > fix(imgParams.nz0 / 2))
    if ((fix(imgParams.nz0 / 2) - 1) > 0)
        zdistcutoff = fix(imgParams.nz0 / 2) - 1;
    else
        zdistcutoff = 0;
    end
end
fprintf("order2=%d, rdistcutoff=%d, zdistcutoff=%f\n", order2, rdistcutoff, zdistcutoff);


kx = k0x * (order2 - order1);
ky = k0y * (order2 - order1);
%{
First time call:
    kx = (100.438431) * (2 - 0) = 200.876862 pixels;	
	ky = (-0.803525 pixels) * (2 - 0) = -1.60705 pixels;
Second time call:
    kx = (100.438431) * (1 - 0) = 100.438431 pixels;
    ky = (-0.803525 pixels) * (1 - 0) = -0.803525 pixels;
%}
otfcutoff = 0.007;
% otfcutoff = 0.6;

if (order1 == 0)
    band1re = bands(:, :, :, 1);
    band1im = 0;
else
    band1re = bands(:, :, :, order1 * 2);
    band1im = bands(:, :, :, order1 * 2 + 1);
end

band2re = bands(:, :, :, order2 * 2);
band2im = bands(:, :, :, order2 * 2 + 1);

% Generate the overlap arrays
overlap0 = makeOverlaps0Kernel(order1, order2, kx, ky, rdistcutoff, zdistcutoff, otfcutoff, order0_2_factor, ...
    band1im, band1re, overlap0, otf, direction, params, imgParams);

overlap1 = makeOverlaps1Kernel(order1, order2, kx, ky, rdistcutoff, zdistcutoff, otfcutoff, order0_2_factor, ...
    band2im, band2re, overlap1, otf, direction, params, imgParams);

% In-place ifft transform of overlap0
overlap0 = ifftn(overlap0);
% In-place ifft transform of overlap1
overlap1 = ifftn(overlap1);
    
end

%% to calculate "m_reconData.savedBands[direction][order = 0] multiplied with order2 otf shifted by (kx, ky)", and saved into overlap0; that is: D0(k)OTF2(k + 2P);
function overlap0 = makeOverlaps0Kernel(order1, order2, kx, ky, rdistcutoff, zdistcutoff, otfcutoff, order0_2_factor,...
    band1im, band1re, overlap0, otf, direction, params, imgParams)

x1 = 0 : (imgParams.nx - 1);
y1 = 0 : (imgParams.ny - 1);
z0 = -zdistcutoff : zdistcutoff;   % z0 = -20:20; if 2D-SIM, z0 = 0 - 0 = 0; Suppose z0 = 10;
x1(x1 >= (imgParams.nx / 2)) = x1(x1 >= (imgParams.nx / 2)) - imgParams.nx;		% x1 = -256:-1:0:255 = -256:255;
y1(y1 >= (imgParams.ny / 2)) = y1(y1 >= (imgParams.ny / 2)) - imgParams.ny;		% y1 = -256:-1:0:255 = -256:255;
% (x1, y1) now is the "relative lateral index to center point (0, 0)" in image after fftshift;
x12 = x1 + kx;
y12 = y1 + ky;
% (x12, y12) now is the "relative lateral index to (-kx, -ky) = (-200.876862, 1.60705)" in image;

otf1 = otf.getOtfMatrix3D(order1, direction, imgParams.nx, imgParams.ny, 0, 0, imgParams.krscale, zdistcutoff, imgParams.kzscale);
otf12 = otf.getOtfMatrix3D(order2, direction, imgParams.nx, imgParams.ny, kx, ky, imgParams.krscale, zdistcutoff, imgParams.kzscale);
% otf1 & otf12 now is a matrix with size [ny, nx, 2 * zdistcutoff + 1]

% rdist1 is the distance from each pixel to the corner point [0, 0][0, 511][511, 0][511, 511]
% or say center point [256, 256] after fftshift;
[X1, Y1, ~] = meshgrid(x1, y1, z0);
rdist1 = sqrt(X1 .* X1 + Y1 .* Y1);

% rdist12 is the distacne from each pixel to [-kx, -ky] = (-200.876862, 1.60705) after fftshift;
[X12, Y12, Z0] = meshgrid(x12, y12, z0);
rdist12 = sqrt(X12 .* X12 + Y12 .* Y12);

overlap_mask = zeros(size(otf1), 'single');
overlap_mask((rdist1 <= rdistcutoff) & (rdist12 <= rdistcutoff) & (abs(otf1) > otfcutoff) & (abs(otf12) * order0_2_factor > otfcutoff)) = 1.0;
% overlap_mask((rdist1 <= (rdistcutoff * otfcutoff)) & (rdist12 <= (rdistcutoff * otfcutoff))) = 1.0;
% if m_myParams.bNoKz0 == true, kz = 0 plane will not be used in modamp fit and assemblerealspace()
if (params.bNoKz0)
    overlap_mask(Z0 == 0) = 0.0;
end
otf1 = otf1 .* overlap_mask;
otf12 = otf12 .* overlap_mask;


z = mod((z0 + imgParams.nz0), imgParams.nz0);     % z0 = [-20:20] --> [80:99]...[0:20] = z
val1re = band1re(:, :, z + 1);
val1im = complex(single(0.0), single(0.0));
if (order1 > 0)
    % Because order1 = fitorder1 = 0, val1im is always 0;
    val1im = band1im(:, :, z + 1);
end


% "root" here is more like a "normalization factor",
% but I suspect that it should be written as: root = abs(otf1) * abs(otf12);
root = sqrt(abs(otf1) .* abs(otf1) + abs(otf12) .* abs(otf12));
% Noise suppresion better?
fact = otf12;
fact(root ~=0) = fact(root ~=0) ./ root(root ~=0);
val1re = val1re .* fact;
if (order1 > 0)
    val1im = val1im .* fact;
end

overlap0(:, :, :) = 0;
overlap0(:, :, z + 1) = val1re - 1j .* val1im;

end

%% to calculate "m_reconData.savedBands[direction][order = -2] multiplied with order0 otf shifted by (-kx, -ky)", and saved into overlap1; that is: D+2(k)OTF0(k - 2P);
%% to calculate "m_reconData.savedBands[direction][order = -1] multiplied with order0 otf shifted by (-kx, -ky)", and saved into overlap1; that is: D+1(k)OTF0(k - P);
function overlap1 = makeOverlaps1Kernel(order1, order2, kx, ky, rdistcutoff, zdistcutoff, otfcutoff, order0_2_factor,...
    band2im, band2re, overlap1, otf, direction, params, imgParams)

x1 = 0 : (imgParams.nx - 1);
y1 = 0 : (imgParams.ny - 1);
z0 = -zdistcutoff : zdistcutoff;   % z0 = -20:20; if 2D-SIM, z0 = 0 - 0 = 0; Suppose z0 = 10;
x1(x1 >= (imgParams.nx / 2)) = x1(x1 >= (imgParams.nx / 2)) - imgParams.nx;		% x1 = -256:-1:0:255 = -256:255;
y1(y1 >= (imgParams.ny / 2)) = y1(y1 >= (imgParams.ny / 2)) - imgParams.ny;		% y1 = -256:-1:0:255 = -256:255;
% (x1, y1) now is the "relative lateral index to center point (0, 0)" in image after fftshift;

x21 = x1 - kx;
y21 = y1 - ky;
% (x21, y21) now is the "relative lateral index to (kx, ky)" in image;

otf2 = otf.getOtfMatrix3D(order2, direction, imgParams.nx, imgParams.ny, 0, 0, imgParams.krscale, zdistcutoff, imgParams.kzscale);
otf21 = otf.getOtfMatrix3D(order1, direction, imgParams.nx, imgParams.ny, -kx, -ky, imgParams.krscale, zdistcutoff, imgParams.kzscale);
% otf1 & otf12 now is a matrix with size [ny, nx, 2 * zdistcutoff + 1]

% rdist1 is the distance from each pixel to the corner point [0, 0][0, 511][511, 0][511, 511]
% or say center point [256, 256] after fftshift;
[X1, Y1, ~] = meshgrid(x1, y1, z0);
rdist1 = sqrt(X1 .* X1 + Y1 .* Y1);

% rdist12 is the distacne from each pixel to [kx, ky] = (200.876862, -1.60705) after fftshift;
[X21, Y21, Z0] = meshgrid(x21, y21, z0);
rdist21 = sqrt(X21 .* X21 + Y21 .* Y21);

overlap_mask = zeros(size(otf2), 'single');
overlap_mask((rdist1 <= rdistcutoff) & (rdist21 <= rdistcutoff) & (abs(otf2) * order0_2_factor > otfcutoff) & (abs(otf21) > otfcutoff)) = 1.0;
% overlap_mask((rdist1 <= (rdistcutoff * otfcutoff)) & (rdist21 <= (rdistcutoff * otfcutoff))) = 1.0;
% if m_myParams.bNoKz0 == true, kz = 0 plane will not be used in modamp fit and assemblerealspace()
if params.bNoKz0
    overlap_mask(Z0 == 0) = 0.0;
end
otf2 = otf2 .* overlap_mask;
otf21 = otf21 .* overlap_mask;


z = mod((z0 + imgParams.nz0), imgParams.nz0);     % z0 = [-20:20] --> [80:99]...[0:20] = z
val2re = band2re(:, :, z + 1);
val2im = band2im(:, :, z + 1);

% "root" here is more like a "normalization factor",
% but I suspect that it should be written as: root = abs(otf1) * abs(otf12);
root = sqrt(abs(otf2) .* abs(otf2) + abs(otf21) .* abs(otf21));
% Noise suppresion better?
fact = otf21;
fact(root ~=0) = fact(root ~=0) ./ root(root ~=0);
val2re = val2re .* fact;
val2im = val2im .* fact;

overlap1(:, :, :) = 0;
overlap1(:, :, z + 1) = val2re - 1j .* val2im;

end

%% to determine the maximum of cross-correlation matrix intensities;
function peak = findpeak(array, sizex, sizey, peak)

% [xcent, ycent] is the index position of largest intensity value;
big = max(array, [], 'all');
[ycent, xcent] = find(array == big);

% a1 is the one-left value of "big";
if (xcent == 1)
    a1 = array(ycent, sizex);
else
    a1 = array(ycent, xcent - 1);
end

% a2 = array[ycent * sizex + xcent];
a2 = big;

% a3 is the one-right value of "big";
if (xcent == sizex)
    a3 = array(ycent, 1);
else
    a3 = array(ycent, xcent + 1);
end

peak.x = fitparabola(a1, a2, a3) + xcent - 1;

% a1 is the one-above value of "big";
if (ycent == 1)
    a1 = array(sizey, xcent);
else
    a1 = array(ycent - 1, xcent);
end

% a2 is the value of "big" itself;
a2 = big;

% a3 is the one-below value of "big";
if (ycent == sizey)
    a3 = array(1, xcent);
else
    a3 = array(ycent + 1, xcent);
end
peak.y = fitparabola(a1, a2, a3) + ycent - 1;

%{
    Since xcent, ycent are within range [1, 512]; fitparabola() return is within range (-1, 1); 
    The return result of peak.x and peak.y should be within (-1, 512);
%}
        
end

%% to localize the sub-pixel position of parabolic fitting maximum;
function peak =  fitparabola(a1, a2, a3)
% INPUT VARIABLES
%   a1:
%   a2:
%   a3:
% OUTPUT VARIABLE
%   peak: x-value of the max (or min) of the parabola
slope = single(0.5 * (a3 - a1));        % the slope at (x=0).
curve = single((a3 + a1) - 2 * a2);     % (a3-a2)-(a2-a1). The change in slope per unit of x.
if (curve == 0)
    disp(['no peak: a1= ', num2str(a1), ', a2=', num2str(a2), ', a3=%f', num2str(a3), ', slope=', num2str(slope), ' curvature=', num2str(curve)]);
    peak = single(0.0);
    return;
end

peak = -slope / curve;          % the x value where slope = 0
if (peak > 1.5 || peak < -1.5)
    fprintf("bad peak position: a1=%f, a2=%f, a3=%f, slope=%f, curvature=%f, peak=%f\n", a1, a2, a3, slope, curve, peak);
    peak = single(0.0);
    return;
end

end

%% Save overlapped bands fft data image
function saveOverlappedbands(params, imgParams, direction, order, overlap0, overlap1)

% Write cross-correlation result image
crosscorr_c = sum(overlap1 .* conj(overlap0), 3);
crosscorr_c = fft2(crosscorr_c);
crosscorr_c = fftshift(abs(crosscorr_c) .^ 2);
fileID = fopen(strcat(params.fileOverlaps, "_dir_", num2str(direction), "_order_", num2str(order), "_xcorr.tif"), 'w+');
[Height, Width] = size(crosscorr_c);
minval = min(crosscorr_c, [], 'all');
maxval = max(crosscorr_c, [], 'all');
header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(crosscorr_c), 1, 1, 1, imgParams.dz, minval, maxval, imgParams.dy);
crosscorr_c = permute(crosscorr_c, [2 1]);
crosscorr_c = reshape(crosscorr_c, [1, Height * Width]);

[~, ~, system_endian] = computer;
if system_endian ~= header.endian
    crosscorr_c = swapbytes(crosscorr_c);
end
crosscorr_c = typecast(crosscorr_c, 'uint8');
fwrite(fileID, crosscorr_c, 'uint8');
fclose(fileID);


% Write overlap0
overlap0 = log(abs(fftshift(fftn(overlap0))));
fileID = fopen(strcat(params.fileOverlaps, "_dir_", num2str(direction), "_order_", num2str(order), "_overlap0.tif"), 'w+');
[Height, Width, Depth] = size(overlap0);
minval = min(overlap0, [], 'all');
maxval = max(overlap0, [], 'all');
header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(overlap0), 1, Depth, 1, imgParams.dz, minval, maxval, imgParams.dy);
overlap0 = permute(overlap0, [2 1 3]);
overlap0 = reshape(overlap0, [1, Height * Width * Depth]);

if system_endian ~= header.endian
    overlap0 = swapbytes(overlap0);
end
overlap0 = typecast(overlap0, 'uint8');
fwrite(fileID, overlap0, 'uint8');
fclose(fileID);

% Write overlap1
overlap1 = log(abs(fftshift(fftn(overlap1))));
fileID = fopen(strcat(params.fileOverlaps, "_dir_", num2str(direction), "_order_", num2str(order), "_overlap1.tif"), 'w+');
[Height, Width, Depth] = size(overlap1);
minval = min(overlap1, [], 'all');
maxval = max(overlap1, [], 'all');
header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(overlap1), 1, Depth, 1, imgParams.dz, minval, maxval, imgParams.dy);
overlap1 = permute(overlap1, [2 1 3]);
overlap1 = reshape(overlap1, [1, Height * Width * Depth]);

if system_endian ~= header.endian
    overlap1 = swapbytes(overlap1);
end
overlap1 = typecast(overlap1, 'uint8');
fwrite(fileID, overlap1, 'uint8');
fclose(fileID);

end

%% (1)to modify / update m_reconData.k0[direction].x and y; (2)to compute m_reconData.amp[direction][order = 1] and [order = 2];
function [k0, amps, overlap0, overlap1] = fitk0andmodamps(bands, overlap0, overlap1, k0, amps, otf, direction, dir_overlap, params, imgParams)
%{
    bands = obj.m_reconData.savedBands(:, :, :, :, direction)
    overlap0 = class: 'single complex'; size:[m_imgParams.ny, m_imgParams.nx, m_imgParams.nz]
    k0 = obj.m_reconData.k0(direction); Suppose direction = 1, k0.x = m_reconData.k0guess(1).x = 100.438431; k0.y = -0.803525
    amps = obj.m_reconData.amp(:, direction)
    otf = obj.otf;
    direction = dir_;
    params = obj.m_myParams;
    imgParams = obj.m_imgParams;
%}

fitorder1 = 0;

if (imgParams.nz0 > 1)
    fitorder2 = 2;    
else
    fitorder2 = 1;
end

k0mag = sqrt(k0.x * k0.x + k0.y * k0.y);
k0angle = atan2(k0.y, k0.x);

% recalculate the overlap arrays at least this first time
redoarrays = (params.recalcarrays >= 1);
x2 = k0angle;
[~, amp2, overlap0, overlap1] = getmodamp(bands, overlap0, overlap1, fitorder1, fitorder2, ...
    k0angle, k0mag, redoarrays, otf, direction, params, imgParams, false);

% recalculate the overlap arrays every time only if recalcarrays >= 3
redoarrays = (params.recalcarrays >= 3);
deltaangle = 0.001;
deltamag = 0.1;
angle = k0angle + deltaangle;
x3 = angle;
[~, amp3, overlap0, overlap1] = getmodamp(bands, overlap0, overlap1, fitorder1, fitorder2, ...
    angle, k0mag, redoarrays, otf, direction, params, imgParams, false);

if (amp3 > amp2)
    while (amp3 > amp2)
        amp1 = amp2;                    % amp1 is the previous modulation factor square;
        x1 = x2;                        % x1 is the previous angle we compute amp2;
        amp2 = amp3;                    % amp2 is the current modulation factor square;
        x2 = x3;                        % x2 is the current angle;
        
        angle = angle + deltaangle;
        x3 = angle;                     % x3 is the new angle we are going to use;
        % Recalculate the modulation factor square amp3 to see whether it is still larger than amp2;
        [~, amp3, overlap0, overlap1] = getmodamp(bands, overlap0, overlap1, fitorder1, fitorder2, ...
            angle, k0mag, redoarrays, otf, direction, params, imgParams, false);
    end
    % After this while loop, amp2 stores the largest modulation factor square; x2 stores the corresponding angle;
else
    % if the modulation factor amp3 <= amp2, go back to k0angle;
    angle = k0angle;
    a = amp3;
    amp3 = amp2;                        % amp3 is the previous modulation factor square (before 0.001 rad adding);
    amp2 = a;                           % amp2 is the current modulation factor square (after 0.001 rad adding);
    a = x3;
    x3 = x2;                            % x3 is the previous angle we compute amp2 (before 0.001 rad increase);
    x2 = a;                             % x2 is the current angle (after 0.001 rad adding);
    while (amp3 > amp2)
        amp1 = amp2;                    % amp1 is the previous modulation factor square;
        x1 = x2;                        % x1 is the previous angle we compute amp2;
        amp2 = amp3;                    % amp2 is the current modulation factor square;
        x2 = x3;                        % x2 is the current angle;
        
        angle = angle - deltaangle;
        x3 = angle;                     % x3 is the new angle we are going to use;
        % Recalculate the modulation factor square amp3 to see whether it is still larger than amp2;
        [~, amp3, overlap0, overlap1] = getmodamp(bands, overlap0, overlap1, fitorder1, fitorder2, ...
            angle, k0mag, redoarrays, otf, direction, params, imgParams, false);
    end
end
% x2 is within [x1:x3]; amp2 is within [amp1:amp3];
angle = fitxyparabola(x1, amp1, x2, amp2, x3, amp3);    % this should be a good angle.


%%%%%   now search for optimum magnitude, at this angle   %%%%%
x2 = k0mag;
[~, amp2, overlap0, overlap1] = getmodamp(bands, overlap0, overlap1, fitorder1, fitorder2, ...
    angle, k0mag, redoarrays, otf, direction, params, imgParams, false);

mag = k0mag + deltamag;
x3 = mag;
[~, amp3, overlap0, overlap1] = getmodamp(bands, overlap0, overlap1, fitorder1, fitorder2, ...
    angle, mag, redoarrays, otf, direction, params, imgParams, false);
if (amp3 > amp2)
    while (amp3 > amp2)
        amp1 = amp2;
        x1 = x2;
        amp2 = amp3;
        x2 = x3;
        
        mag = mag + deltamag;
        x3 = mag;
        [~, amp3, overlap0, overlap1] = getmodamp(bands, overlap0, overlap1, fitorder1, fitorder2, ...
            angle, mag, redoarrays, otf, direction, params, imgParams, false);
    end
else
    mag = k0mag;
    a = amp3;
    amp3 = amp2;
    amp2 = a;
    a = x3;
    x3 = x2;
    x2 = a;
    while (amp3 > amp2)
        amp1 = amp2;
        x1 = x2;
        amp2 = amp3;
        x2 = x3;
        
        mag = mag - deltamag;
        x3 = mag;
        [~, amp3, overlap0, overlap1] = getmodamp(bands, overlap0, overlap1, fitorder1, fitorder2, ...
            angle, mag, redoarrays, otf, direction, params, imgParams, false);
    end    
end
% x2 is within [x1:x3]; amp2 is within [amp1:amp3];
mag = fitxyparabola(x1, amp1, x2, amp2, x3, amp3);  % this should be a good magnitude.
% if we were perfectionist we'd iterate for angle again now

fprintf("\n");
fprintf("Optimum modulation amplitude:\n");
redoarrays = (params.recalcarrays >= 2);
[modamp, ~, overlap0, overlap1] = getmodamp(bands, overlap0, overlap1, fitorder1, fitorder2, ...
    angle, mag, redoarrays, otf, direction, params, imgParams, true);
% one last time, to find the modamp at the optimum k0

dk = (1 / (imgParams.ny * imgParams.dy));   % inverse microns per pixel in data
fprintf("Optimum k0 angle=%f, length=%f, spacing=%f microns\n\n", angle, mag, 1.0 / (mag * dk));

if (params.bSaveOverlaps)
    fprintf("Saving overlapped bands\n");
    saveOverlappedbands(params, imgParams, dir_overlap, fitorder2, overlap0, overlap1);
    fprintf("Overlapped bands saved\n\n");
end

k0.x = mag * cos(angle);
k0.y = mag * sin(angle);
amps(fitorder2 + 1) = modamp;

% finally find the modamp for the other orders
redoarrays = true;
if (imgParams.nz0 == 1)
    for order = 2 : (params.norders - 1)
        % assuming that "angle" and "mag" remain the same for every adjacent pair of bands within one direction
        [modamp, ~, overlap0, overlap1] = getmodamp(bands, overlap0, overlap1, order - 1, order, ...
            angle, mag, redoarrays, otf, direction, params, imgParams, true);
        amps(order + 1) = modamp;
        
        fprintf("\n");
        if (params.bSaveOverlaps)
            fprintf("Saving overlapped bands\n");
            saveOverlappedbands(params, imgParams, dir_overlap, order, overlap0, overlap1);
            fprintf("Overlapped bands saved\n");
        end
    end
else
    for order = 1 : (params.norders - 1)
        if (order ~= fitorder2)
            [modamp, ~, overlap0, overlap1] = getmodamp(bands, overlap0, overlap1, 0, order, ...
                angle, mag, redoarrays, otf, direction, params, imgParams, true);
            amps(order + 1) = modamp;
            
            fprintf("\n");
            if (params.bSaveOverlaps)
                fprintf("Saving overlapped bands\n");
                saveOverlappedbands(params, imgParams, dir_overlap, order, overlap0, overlap1);
                fprintf("Overlapped bands saved\n");
            end
        end
    end
end
fprintf("\n");

end

%% to compute modulation factor square: amp2 = modamp.x * modamp.x + modamp.y * modamp.y;
function [modamp, amp2, overlap0, overlap1] = getmodamp(bands, overlap0, overlap1, order1, order2, kangle, klength, redoarrays, otf, direction, params, imgParams, bShowDetail)
%{
    order1 = 0;
	order2 = 2;
	or order2 = 1 to compute order = 1 modulation factor;
    kangle = k0angle = atan2(k0.y, k0.x);
	klength = sqrt(k0.x * k0.x + k0.y * k0.y);
%}
k1 = struct('x', single(0.0), 'y', single(0.0));
k1.x = klength * cos(kangle);
k1.y = klength * sin(kangle);

[modamp, amp_inv, amp_combo, corr_coef, overlap0, overlap1] = findrealspacemodamp(bands, overlap0, overlap1, order1, order2, k1, redoarrays, ...
    otf, direction, params, imgParams);

amp2 = abs(modamp)^2;
fprintf(" In getmodamp: angle=%f, mag=%f, amp=%f, phase=%f\n", kangle, klength, sqrt(amp2), angle(modamp));

if (bShowDetail)
    fprintf(" Reverse modamp is: amp=%f, phase=%f\n", 1.0 / abs(amp_inv), -angle(amp_inv));
    fprintf(" Combined modamp is: amp=%f, phase=%f\n", abs(amp_combo), angle(amp_combo));
    fprintf(" Correlation coefficient is: %f\n", corr_coef);	
end

end

%% to compute modamp, amp_inv & amp_combo;
function [modamp1, modamp2, modamp3, corr_coef, overlap0, overlap1] = findrealspacemodamp(bands, overlap0, overlap1, order1, order2, k0, redoarrays, otf, direction, params, imgParams)
%{
    order1 = 0;
	order2 = 2;
	or order2 = 1 to compute order = 1 modulation factor;
    redoarrays = true;
%}

if (redoarrays)
    % make arrays that contain only the overlapping parts of fourier space. Otf-equalize there, set to zero elsewhere
    [overlap0, overlap1] = makeoverlaps(bands, overlap0, overlap1, order1, order2, k0.x, k0.y, otf, direction, params, imgParams);
    %{
        (1) calculate overlap0 = "m_reconData.savedBands[direction][order = 0] multiplied with order2 otf shifted by (k0.x, k0.y)";		that is: D0(k)OTF2(k + 2P)
		(2) calculate overlap1 = "m_reconData.savedBands[direction][order = +2] multiplied with order0 otf shifted by (-k0.x, -k0.y)";	that is: D+2(k)OTF0(k - 2P);
		(3) apply in-place ifft transform of overlap0 and overlap1;
        Noted that: (k0.x, k0.y) has already been the initial estimate of modulation wave vector m_reconData.k0[direction] by cross-correlation here; 
    %}
end

kx = k0.x * (order2 - order1);
ky = k0.y * (order2 - order1);

i = 0 : (imgParams.nx - 1);
j = 0 : (imgParams.ny - 1);
[X, Y] = meshgrid(i, j);

% phase_shift here is the phase of illumination pattern (kx, ky)
phase_shift = single(2.0 * pi * (X .* kx ./ imgParams.nx + Y .* ky ./ imgParams.ny));

% overlap1 * exp(-1j * phase_shift) = ifft transform of "m_reconData.savedBands[direction][order = +2] shifted by (k0.x, k0.y) multiplied with order0 otf"; that is: ifft[D+2(k + 2P)OTF0(k)];
XStarYFR = sum(sum((conj(overlap0) .* overlap1), 3) .* exp(-1j .* phase_shift), 'all');
sumXMagFR = sum(abs(overlap0) .^ 2, 'all');
sumYMagFR = sum(abs(overlap1) .^ 2, 'all');
%{
Now, XStarYFR is the sum of: conj(overlap0) * overlap1 * exp(-j * phase_shift);
	//		sumXMagFR is the sum of: overlap0 * conj(overlap0);
	//		sumYMagFR is the sum of: overlap1 * conj(overlap1)
%}

%{
    Parseval's theorem --->		Sum(f(x) * conj(g(x))) = Sum(fft[f(x)] * conj(fft[g(x)]));
	Plancherel's theorem --->	Sum(f(x) * conj(f(x))) = Sum(fft[f(x)] * conj(fft[f(x)]));
	In our application, we want to calculate the sum of conjugate of overlap0? "m_reconData.savedBands[direction][order = 0] multiplied with order2 otf shifted by (-kx, -ky)"
		multiplied with overlap1 shifted by (kx, ky): "m_reconData.savedBands[direction][order = minus 2] shifted by (kx, ky) multiplied with order0 otf",
		which are all in Fourier domain;
	So all the Fourier transform should be modified by inverse Fourier transform:		
		Sum(conj(overlap0) * overlap1(shift)) = Sum(conj(ifft[overlap0]) * ifft[overlap1(shift)]) = Sum(conj(ifft[overlap0]) * ifft[overlap1] * exp(j * shift));
	
	That is the reason why we apply in-place ifft transform of overlap0 and overlap1 here.
	After in-place ifft transform of overlap0 and overlap1, Sum = Sum(conj(overlap0) * overlap1 * exp(j * phase_shift)) = XStarYFR;
	
	Similarly, we want to calculated Sum(overlap0 * conj(overlap0)) = Sum(ifft[overlap0] * conj(ifft[overlap0]));
	After in-place ifft transform of overlap0,	Sum = Sum(overlap0 * conj(overlap0));
%}

modamp1 = XStarYFR / sumXMagFR; % modamp
modamp2 = XStarYFR / sumYMagFR; % amp_inv

tan2beta = 2.0 * abs(XStarYFR) / (sumXMagFR - sumYMagFR);
beta = 0.5 * atan(tan2beta);
if (beta < 0.0)
    beta = beta + 0.5 * pi;      % Convert (-pi / 4, pi / 4) to (0, pi / 2);
end

modamp3 = tan(beta) * exp(1j * angle(XStarYFR));
corr_coef = abs(XStarYFR) / sqrt(sumXMagFR * sumYMagFR);    

end

%% to localize the sub-pixel position of parabolic fitting maximum of angle and magnitude;
function peak = fitxyparabola(x1, y1, x2, y2, x3, y3)
% x2 is within [x1:x3]; y2 is within [y1:y3];

if (x1 == x2 || x2 == x3 || x3 == x1)
    fprintf("Fit fails; two points are equal: x1=%f, x2=%f, x3=%f\n", x1, x2, x3);
    peak = 0.0;
    return;
end

xbar1 = 0.5 * (x1 + x2);            % middle of x1 and x2
xbar2 = 0.5 * (x2 + x3);            % middle of x2 and x3
slope1 = (y2 - y1) / (x2 - x1);     % the slope at (x = xbar1).
slope2 = (y3 - y2) / (x3 - x2);     % the slope at (x = xbar2).
curve = (slope2 - slope1) / (xbar2 - xbar1);    % The change in slope per unit of x.
if (curve == 0)
    fprintf("Fit fails; no curvature: r1=(%f,%f), r2=(%f,%f), r3=(%f,%f) slope1=%f, slope2=%f, curvature=%f\n", x1, y1, x2, y2, x3, y3, slope1, slope2, curve);
    peak = 0.0;
    return;
end

peak = xbar2 - slope2 / curve;

end





%% Calculate reconstruction parameters, including rdistcutoff, zdistcutoff, apocutoff and zapocutoff
function [params, imgParams] = calculateReconstructionParameters(k0, params, imgParams)

SPOTRATIO  = 0;

% k0 magnitude (for 1st order) in pixels
% The reason why we only use k0 of 1st direction is based on the fact that:
% all 3 directions should not have too much variance
k0pix = sqrt(k0(1).x ^ 2 + k0(1).y ^ 2);
k0pix = k0pix * (params.norders - 1);
% k0 magnitude (for highest order) in inverse microns
k0mag = k0pix * imgParams.dkr;


% emission wavelength in the sample, in microns
lambdaem = (imgParams.wave(1) / params.nimm) / 1.0e3;
% 0.93 approximates a typical lambdaexc/lambdaem
% lambdaexc = 0.93 * lambdaem;
lambdaexc = (params.lambda / params.nimm) / 1.0e3;
% lambdaexc_medium is to determine standing wave position in frequency domain
lambdaexc_medium = params.lambda / params.nmed / 1.0e3;
% alpha is the aperture angle of objectives
alpha = asin(params.na / params.nimm);
% beta is the angle of center of side illumination beams
beta = asin(k0mag / (2.0 / lambdaexc_medium));
% betamin is the angle of inner edge of side illumination beams (for multi-mode laser)
betamin = asin(k0mag / (2.0 / lambdaexc_medium) - sin(alpha) * SPOTRATIO);

% OTF support radial limit in data pixels
imgParams.rdistcutoff = fix((params.na * 2.0 / (imgParams.wave(1) / 1.0e3)) / imgParams.dkr);
if (imgParams.rdistcutoff > fix(imgParams.nx / 2))
    imgParams.rdistcutoff = fix(imgParams.nx / 2);              % rdistcutoff = 210 < 256 = nx / 2;
end

imgParams.zdistcutoff = zeros([params.norders, 1], 'single');

if (~params.bTwolens)
    % OTF support axial limit in data pixels
    imgParams.zdistcutoff(1) = ceil(((1 - cos(alpha)) / lambdaem) / imgParams.dkz);
    
    % approx max axial support limit of the OTF of the high frequency side band
    imgParams.zdistcutoff(params.norders) = fix(1.02 * imgParams.zdistcutoff(1));
    if (params.norders >= 3)
        for order = 1 : (params.norders - 2)
            % axial support limit of the OTF of the medium frequency side bands
            %             imgParams.zdistcutoff(order + 1) = fix((1 + lambdaem / lambdaexc) * imgParams.zdistcutoff(1));
            imgParams.zdistcutoff(order + 1) = fix(((1 - cos(alpha)) / lambdaem + (1 - cos(betamin)) / lambdaexc) / imgParams.dkz);
        end
    end
else
    % 1.02 is just a safety margin
    imgParams.zdistcutoff(1) = ceil(1.02 * ((1 - cos(alpha)) / lambdaem + 2 / lambdaexc_medium) / imgParams.dkz);
    % approx max axial support limit of the OTF of the high frequency side band
    imgParams.zdistcutoff(params.norders) = fix(1.02 * ceil(((1 - cos(alpha)) / lambdaem) / imgParams.dkz));
    
    if (params.norders == 3)
        % axial support limit of the OTF of the medium frequency side band
        imgParams.zdistcutoff(2) = ceil(1.02 * ((1 - cos(alpha)) / lambdaem + (1 + cos(betamin)) / lambdaexc) / imgParams.dkz);
    elseif (params.norders > 3)
        for order = 1 : (params.norders - 1)
            a = order / (params.norders - 1);
            % axial support limit of the OTF of the medium frequency side bands, 1.1 is a blur margin */
            imgParams.zdistcutoff(order + 1) = fix(1.1 * ((1 - a) * imgParams.zdistcutoff(1) + a * imgParams.zdistcutoff(params.norders)));
        end
    end
end

for order = 1 : params.norders
    if (imgParams.zdistcutoff(order) >= fix(imgParams.nz0 / 2))
        if ((fix(imgParams.nz0 / 2) - 1) > 0)
            imgParams.zdistcutoff(order) = fix(imgParams.nz0 / 2) - 1;    % 3D-SIM
        else
            imgParams.zdistcutoff(order) = 0;                             % 2D-SIM
        end
    end
end


% apocutoff is the OTF support radial limit after reconstruction
% imgParams.apocutoff = imgParams.rdistcutoff + k0pix * (params.norders - 1);
imgParams.apocutoff = imgParams.rdistcutoff + k0pix;

if (params.bTwolens)
    imgParams.zapocutoff = imgParams.zdistcutoff(1);        % If I5S, OTF support axial limit after reconstruction takes from 0th order
else
    imgParams.zapocutoff = imgParams.zdistcutoff(2);        % If 3D-SIM, OTF support axial limit after reconstruction takes from 1st order
end

end



%% to apply RL deconvolution filter
function [fullResult, otfSim] = filterbands_RL(dir, bands, k0, amp, noiseVarFactors, otf, params, imgParams)

% Get the "genuine" seprated bands
separate = zeros(size(bands), 'single');
separate(:, :, :, 1) = bands(:, :, :, 1);
for order = 1 : (params.norders - 1)
    % Minus bands
    separate(:, :, :, 2 * order)     = (bands(:, :, :, 2 * order) + 1j * bands(:, :, :, 2 * order + 1)) ./ conj(amp(order + 1, dir));
    % Plus bands
    separate(:, :, :, 2 * order + 1) = (bands(:, :, :, 2 * order) - 1j * bands(:, :, :, 2 * order + 1)) ./ amp(order + 1, dir);
end

dir_otf = 1;
if (params.bOneOTFperAngle)
    dir_otf = dir;
end

% Shift separated bands and OTF to correct position
shifted = zeros([params.zoomfact * imgParams.ny, params.zoomfact * imgParams.nx, params.z_zoom * imgParams.nz0, params.nphases], 'single');
otfSim = zeros([params.zoomfact * imgParams.ny, params.zoomfact * imgParams.nx, params.z_zoom * imgParams.nz0], 'single');

for order = -(params.norders - 1) : (params.norders - 1)
    kx = order * k0(dir).x;
    ky = order * k0(dir).y;
    fprintf("dir %d: shift order %d to (%.3f, %.3f)\n", dir, order, kx, ky);
    
    % Shift separated bands
    if (order == 0)
        % band 0 is DC, so does not need shifting, only a bigger matrix
        shifted(:, :, :, 1) = placeFreq(separate(:, :, :, 1), shifted(:, :, :, 1));
    elseif (order < 0)
        % higher bands need shifting
        % first copy to larger matrix, then Fourier shift
        shifted(:, :, :, 2 * abs(order)) = placeFreq(separate(:, :, :, 2 * abs(order)), shifted(:, :, :, 2 * abs(order)));
        shifted(:, :, :, 2 * abs(order)) = fourierShift(shifted(:, :, :, 2 * abs(order)), kx, ky);
    elseif (order > 0)
        shifted(:, :, :, 2 * abs(order) + 1) = placeFreq(separate(:, :, :, 2 * abs(order) + 1), shifted(:, :, :, 2 * abs(order) + 1));
        shifted(:, :, :, 2 * abs(order) + 1) = fourierShift(shifted(:, :, :, 2 * abs(order) + 1), kx, ky);
    end
    
    % Shift OTF
    OTF = zeros(size(shifted(:, :, :, 2 * abs(order) + 1)), 'single');
    otfval = otf.getOtfMatrix3D(abs(order), dir_otf, params.zoomfact * imgParams.nx, params.zoomfact * imgParams.ny, kx, ky, imgParams.krscale, imgParams.zdistcutoff(abs(order) + 1), imgParams.kzscale);
%     if (order == 0)
%         otfval = otfval / noiseVarFactors(dir, abs(order) + 1);
%     else
%         dampen_matrix = otf.getAttenuationMatrix3D(abs(order), params.zoomfact * imgParams.nx, params.zoomfact * imgParams.ny, kx, ky, imgParams.zdistcutoff(abs(order) + 1), imgParams.rdistcutoff, params);
%         otfval = otfval .* dampen_matrix / noiseVarFactors(dir, abs(order) + 1);
%     end
    
    %     otfval = otfval ./ max(abs(otfval), [], 'all');
    
    z0 = -imgParams.zdistcutoff(abs(order) + 1) : imgParams.zdistcutoff(abs(order) + 1);
    iz = mod((z0 + params.z_zoom * imgParams.nz0), params.z_zoom * imgParams.nz0);
    OTF(:, :, iz + 1) = otfval;
    
    if (order == 0)
        OTF = OTF / params.ndirs;
    end
    otfSim = otfSim + OTF;    
end

fprintf("\n");

shifted(:, :, :, 1) = shifted(:, :, :, 1) / params.ndirs;
fullResult = sum(shifted, 4);

end

%% to pre-filter separated bands for one directions;
function filterbands_single_direction(it, dir, bands, k0, amp, noiseVarFactors, otf, params, imgParams)
%{
    dir = 1, 2, 3; Suppose dir = 1;
    bands = obj.m_reconData.savedBands(:, :, :, :, direction);
    k0 =  obj.m_reconData.k0;
    amp = obj.m_reconData.amp;
    noiseVarFactors = obj.m_reconData.noiseVarFactors;
    otf = obj.otf;
    params = obj.m_myParams;
    imgParams = obj.m_imgParams;
%}
% Option #1: multiply separated bands with conjugated OTF and phase factors, both unshifted

fact = params.explodefact / 0.5;

dir_otf = 1;
if (params.bOneOTFperAngle)
    dir_otf = dir;
end

bigbuffer = zeros([params.zoomfact * imgParams.ny, params.zoomfact * imgParams.nx, params.z_zoom * imgParams.nz0], 'single');
outbuffer = zeros([params.zoomfact * imgParams.ny, params.zoomfact * imgParams.nx, params.z_zoom * imgParams.nz0], 'single');

for order = 0 : (params.norders - 1)
    fprintf("dir %d: filtering order %d\n", dir, order);
    kx = order * k0(dir).x;
    ky = order * k0(dir).y;
    
    
    % otf1 is the order = 0, 1, 2 OTF value corresponding to the current (x1, y1, z0) or say "relative index to center point (0, 0, 0)" in image;
    otf1 = otf.getOtfMatrix3D(order, dir_otf, imgParams.nx, imgParams.ny, 0, 0, imgParams.krscale, imgParams.zdistcutoff(order + 1), imgParams.kzscale);
    % otf1: OTF(order = m)(k)
    weight = abs(otf1) .^ 2;
    if (order ~= 0)
        weight = weight .* abs(amp(order + 1, dir)) ^ 2;
    end
    
    dampen_matrix = otf.getAttenuationMatrix3D(order, imgParams.nx, imgParams.ny, 0, 0, imgParams.zdistcutoff(order + 1), imgParams.rdistcutoff, params);
    dampfact = dampen_matrix / noiseVarFactors(dir, order + 1);
    
    weight = weight .* dampfact;
    sumweight = weight;
    
    % Calculate all-direction Wiener denominator
    for order2 = -(params.norders - 1) : (params.norders - 1)
        % order2 = [-2, -1, 0, 1, 2];
        if (order2 == order)
            continue;
        end
        
        %  bFilteroverlaps is always true except when (during debug) generating an unfiltered exploded view
        if (~params.bFilteroverlaps && ~(order2 == 0 && order == 0))
            continue;
        end
        
        % shifted center of band 2
        kx2 = order2 * k0(dir).x;
        ky2 = order2 * k0(dir).y;
        
        dir_otf = 1;
        if (params.bOneOTFperAngle)
            dir_otf = dir;
        end
        
        % otf2 is the abs(order2) = 0, 1, 2 OTF value corresponding to the current (x2, y2, z0) or say "relative index to center point (0, 0, 0)" in image;
        otf2 = otf.getOtfMatrix3D(abs(order2), dir_otf, imgParams.nx, imgParams.ny, kx2 - kx, ky2 - ky, imgParams.krscale, imgParams.zdistcutoff(order + 1), imgParams.kzscale);
        weight = abs(otf2) .^ 2;
        weight = weight ./ noiseVarFactors(dir, abs(order2) + 1);
        
        amp2mag2 = abs(amp(abs(order2) + 1, dir)) ^ 2;
        if (order2 ~= 0)
            weight = weight .* amp2mag2;
        end
        
        dampen_matrix = otf.getAttenuationMatrix3D(abs(order2), imgParams.nx, imgParams.ny, kx2 - kx, ky2 - ky, imgParams.zdistcutoff(order + 1), imgParams.rdistcutoff, params);
        weight = weight .* dampen_matrix;
        sumweight = sumweight + weight;
    end
    % After the for loop: sumweight = dampfact * (otf1.x ^ 2 + otf1.y ^ 2) * Magnitude[m_reconData.amp[dir][order]] + Sum[dev_suppress(rdist2) * (otf2.x ^ 2 + otf2.y ^ 2) / 1.0f * Magnitude[m_reconData.amp[dir][order2 = 0,1,2]]];
    % 5 orders in total;
    
    sumweight = sumweight + params.wiener ^ 2;
    % Above is the denominator part of the shifted wiener filter function;
    
    scale = dampfact .* conj(otf1);
    % scale = dampfact .* conj(OTF(dir, order) = dampfact .* O(*)(m)(k)
    scale = scale ./ sumweight;
    % Based on 2008 Biophysical Journal paper's defination:
    % scale = conj[OTF(dir, order)] / (sum[abs(OTF(dir, order2) shifted by [-order2 * p(dir) + order * k0(dir)]) .^ 2] + wiener ^ 2);

    
    % Apodization function: A(k - mp(d)) = A(k - order * k0), which is function shifted by -mp(d)
    if (params.apodizeoutput)
        apofact = otf.writeApoFactor(imgParams.nx, imgParams.ny, -kx, -ky, imgParams.zdistcutoff(order + 1), imgParams.apocutoff, imgParams.zapocutoff, params);
        scale = scale .* apofact;
        % Based on 2008 Biophysical Journal paper's defination:
        % scale = conj[OTF(order)] * A(k - order * p(dir)) /
        % (sum[abs(OTF(dir, order2) shifted by [order2 * p(dir) - order * k0(dir)]) .^ 2] + wiener ^ 2);
        % Simply speaking, scale = A(k - order * p(dir)) * conjugate(otf1)/sumweight
    end
    
    z0 = -imgParams.zdistcutoff(order + 1) : imgParams.zdistcutoff(order + 1);  % z0 = -20:0:20
    iz = mod((z0 + imgParams.nz0), imgParams.nz0);          % z0 = [-20:20] --> [80:99]...[0:20] = iz
    
    %{
    separate (for this pixel) the even and odd "bands" into the true plus and minus bands
	apply the scale to the plus band only (this will apply it to the minus band as well by symmetry(?)
    reassemble into the even and odd "bands"
    %}
    band = zeros([imgParams.ny, imgParams.nx, imgParams.nz0], 'single');
    if order == 0
        band(:, :, iz + 1) = bands(:, :, iz + 1, 1) .* scale;
        
        % 0th order
        fprintf("moving centerband\n");
        bigbuffer(:) = 0;
        bigbuffer = placeFreq(band(:, :, :), bigbuffer);
        
        % transform it into real space
        fprintf("re-transforming central band\n");
        bigbuffer = ifftn(bigbuffer);
        
        fprintf("inserting central band\n");
        outbuffer = write_outbuffer_kernel1(bigbuffer, outbuffer);
        
        fprintf("central band assembly completed\n");
    else
        scale = scale .* conj(amp(order + 1, dir));
        % Based on 2008 Biophysical Journal paper's defination:
        % scale = A(k - order * p(dir)) * [c(m)* exp(-1j * phi(m))] * conjugate(otf1)/sumweight
        % conjamp(order + 1) = c(m)* exp(1j * phi(m))        
        band(:, :, iz + 1) = (bands(:, :, iz + 1, 2 * order) - 1j * bands(:, :, iz + 1, 2 * order + 1)) .* scale;
        
        fprintf("moving order %d\n", order);
        bigbuffer(:) = 0;
        bigbuffer = placeFreq(band(:, :, :), bigbuffer);
        
        % transform it into real space
        fprintf("re-transforming side bands\n");
        bigbuffer = ifftn(bigbuffer);
        
        % For 3D, prepare 2D array of sines and cosines first, then loop over z.
        k0x = k0(dir).x * order;
        k0y = k0(dir).y * order;
        x = 0 : (params.zoomfact * imgParams.nx - 1);
        y = 0 : (params.zoomfact * imgParams.ny - 1);
        z = 0 : (params.z_zoom * imgParams.nz0 - 1);
        [X, Y, ~] = meshgrid(x, y, z);
        angle = fact * pi * (X .* k0x / (params.zoomfact * imgParams.nx) + Y .* k0y / (params.zoomfact * imgParams.ny));
        coslookup = cos(-angle);
        sinlookup = sin(-angle);
        
        fprintf("inserting side bands\n");
        outbuffer = write_outbuffer_kernel2(coslookup, sinlookup, bigbuffer, outbuffer);
        
        fprintf("order %d side band assembly completed\n", order);
    end
    fprintf("\n");  

end

saveWienerResult_single_direction(dir, it, outbuffer, params, imgParams);

end

%% to pre-filter separated bands for one directions (shifted);
function fullResult = filterbands_single_direction_shifted(it, dir, bands, k0, amp, noiseVarFactors, otf, params, imgParams)
%{
    dir = 1, 2, 3; Suppose dir = 1;
    bands = obj.m_reconData.savedBands(:, :, :, :, direction);
    k0 =  obj.m_reconData.k0;
    amp = obj.m_reconData.amp;
    noiseVarFactors = obj.m_reconData.noiseVarFactors;
    otf = obj.otf;
    params = obj.m_myParams;
    imgParams = obj.m_imgParams;
%}
% Option #2: multiply separated bands (shifted) with conjugated OTF (shifted) and phase factors

% Get the "genuine" seprated bands
separate = zeros(size(bands), 'single');
separate(:, :, :, 1) = bands(:, :, :, 1);
for order = 1 : (params.norders - 1)
    % Minus bands
    separate(:, :, :, 2 * order)     = bands(:, :, :, 2 * order) + 1j * bands(:, :, :, 2 * order + 1);
    % Plus bands
    separate(:, :, :, 2 * order + 1) = bands(:, :, :, 2 * order) - 1j * bands(:, :, :, 2 * order + 1);
end

dir_otf = 1;
if (params.bOneOTFperAngle)
    dir_otf = dir;
end

% Shift separated bands to correct position
shifted = zeros([params.zoomfact * imgParams.ny, params.zoomfact * imgParams.nx, params.z_zoom * imgParams.nz0, params.nphases], 'single');

for order = -(params.norders - 1) : (params.norders - 1)
    kx = order * k0(dir).x;
    ky = order * k0(dir).y;
    fprintf("dir %d: shift order %d to (%.3f, %.3f)\n", dir, order, kx, ky);
    
    if (order == 0)
        % band 0 is DC, so does not need shifting, only a bigger matrix
        shifted(:, :, :, 1) = placeFreq(separate(:, :, :, 1), shifted(:, :, :, 1));
    elseif (order < 0)
        % higher bands need shifting
        % first copy to larger matrix, then Fourier shift
        shifted(:, :, :, 2 * abs(order)) = placeFreq(separate(:, :, :, 2 * abs(order)), shifted(:, :, :, 2 * abs(order)));
        shifted(:, :, :, 2 * abs(order)) = fourierShift(shifted(:, :, :, 2 * abs(order)), kx, ky);
    elseif (order > 0)
        shifted(:, :, :, 2 * abs(order) + 1) = placeFreq(separate(:, :, :, 2 * abs(order) + 1), shifted(:, :, :, 2 * abs(order) + 1));
        shifted(:, :, :, 2 * abs(order) + 1) = fourierShift(shifted(:, :, :, 2 * abs(order) + 1), kx, ky);
    end
end
fprintf("\n");


for order = -(params.norders - 1) : (params.norders - 1)
    kx = order * k0(dir).x;
    ky = order * k0(dir).y;
    fprintf("dir %d: multiply order %d with conjugated OTF\n", dir, order);
    
    OTF = zeros(size(shifted(:, :, :, 2 * abs(order) + 1)), 'single');
    otfval = otf.getOtfMatrix3D(abs(order), dir_otf, params.zoomfact * imgParams.nx, params.zoomfact * imgParams.ny, kx, ky, imgParams.krscale, imgParams.zdistcutoff(abs(order) + 1), imgParams.kzscale);
    dampen_matrix = otf.getAttenuationMatrix3D(abs(order), params.zoomfact * imgParams.nx, params.zoomfact * imgParams.ny, kx, ky, imgParams.zdistcutoff(abs(order) + 1), imgParams.rdistcutoff, params);
    otfval = otfval .* dampen_matrix / noiseVarFactors(dir, abs(order) + 1);
    
    z0 = -imgParams.zdistcutoff(abs(order) + 1) : imgParams.zdistcutoff(abs(order) + 1);
    iz = mod((z0 + params.z_zoom * imgParams.nz0), params.z_zoom * imgParams.nz0);
    OTF(:, :, iz + 1) = conj(otfval);
    
    if (order == 0)
        shifted(:, :, :, 1) = shifted(:, :, :, 1) .* OTF;
    elseif (order < 0)
        shifted(:, :, :, 2 * abs(order)) = shifted(:, :, :, 2 * abs(order)) .* OTF .* amp(abs(order) + 1, dir);
    elseif (order > 0)
        shifted(:, :, :, 2 * abs(order) + 1) = shifted(:, :, :, 2 * abs(order) + 1) .* OTF .* conj(amp(abs(order) + 1, dir));
    end
end
fprintf("\n");


result = sum(shifted, 4);
fullResult = result;
fDenom = params.wiener ^ 2;
% Calculate per-direction Wiener denominator
for order = -(params.norders - 1) : (params.norders - 1)
    kx = order * k0(dir).x;
    ky = order * k0(dir).y;
    
    OTF = zeros(size(shifted(:, :, :, 1)), 'single');
    otfval = otf.getOtfMatrix3D(abs(order), dir_otf, params.zoomfact * imgParams.nx, params.zoomfact * imgParams.ny, kx, ky, imgParams.krscale, imgParams.zdistcutoff(abs(order) + 1), imgParams.kzscale);
    otfval = abs(otfval) .^2 * abs(amp(abs(order) + 1, dir)) ^ 2 ./ noiseVarFactors(dir, abs(order) + 1);
    dampen_matrix = otf.getAttenuationMatrix3D(abs(order), params.zoomfact * imgParams.nx, params.zoomfact * imgParams.ny, kx, ky, imgParams.zdistcutoff(abs(order) + 1), imgParams.rdistcutoff, params);
    otfval = otfval .* dampen_matrix;
    
    z0 = -imgParams.zdistcutoff(abs(order) + 1) : imgParams.zdistcutoff(abs(order) + 1);
    iz = mod((z0 + params.z_zoom * imgParams.nz0), params.z_zoom * imgParams.nz0);
    OTF(:, :, iz + 1) = otfval;
    fDenom = fDenom + OTF;
end

if (params.apodizeoutput)
    apo = zeros(size(OTF), 'single');
    apofact = otf.writeApoFactor(params.zoomfact * imgParams.nx, params.zoomfact * imgParams.ny, 0, 0, imgParams.zapocutoff, imgParams.apocutoff, imgParams.zapocutoff, params);
    
    z0 = -imgParams.zapocutoff : imgParams.zapocutoff;
    iz = mod((z0 + params.z_zoom * imgParams.nz0), params.z_zoom * imgParams.nz0);
    apo(:, :, iz + 1) = apofact;
    
    result = result .* apo;
end

result = result ./ fDenom;
saveWienerResult_single_direction_shifted(dir, it, result, params, imgParams);

end

%% Place a 3d freq space matrix into a matrix of double size.
function outarray = placeFreq(inarray, outarray)

[ny, nx, nz] = size(inarray);
[ydim, xdim, zdim] = size(outarray);

x = 0 : (nx - 1);
y = 0 : (ny - 1);
z = 0 : (nz - 1);

x(x >= (nx / 2)) = x(x >= (nx / 2)) - nx;       % x = [0:255]...[-256:-1]
y(y >= (ny / 2)) = y(y >= (nx / 2)) - ny;       % y = [0:255]...[-256:-1]
if (nz > 1)
    z(z >= (nz / 2)) = z(z >= (nz / 2)) - nz;   % z = [0:49]...[-50:-1]
end

xout = x;                                   % xout = [0:255]...[-256:-1]
xout(xout < 0) = xout(xout < 0) + xdim;     % xout = [0:255]...[768:1023]
yout = y;                                   % yout = [0:255]...[-256:-1]
yout(yout < 0) = yout(yout < 0) + ydim;     % yout = [0:255]...[768:1023]
zout = z;                                   % zout = [0:49]...[-50:-1]
zout(zout < 0) = zout(zout < 0) + zdim;     % zout = [0:49]...[50:99]

outarray(yout + 1, xout + 1, zout + 1) = inarray;

end

%% Moves freq space data to kx, ky with subpixel precision, by phase multiplication in real space. Basically, just ifft, multiply phases, fft...
function outarray = fourierShift(inarray, kx, ky)

[ny, nx, nz] = size(inarray);
i = 0 : (nx - 1);
j = 0 : (ny - 1);
k = 0 : (nz - 1);
[X, Y, ~] = meshgrid(i, j, k);

% phase_shift here is the phase of illumination pattern (kx, ky)
phase_shift = single(2.0 * pi * (X .* kx ./ nx + Y .* ky ./ ny));

inarray = ifftn(inarray);
outarray = fftn(inarray .* exp(-1j .* phase_shift));

end

%% Save per-direction Wiener filtered result
function saveWienerResult_single_direction(dir, it, image, params, imgParams)

% Wrtie per-direction Wiener filtered fft data image
if (params.bsavefft)
    image_fft = image;
    if ~(params.bfftkeepnegative)
        image_fft(image_fft < 0) = 0;
    end
    image_fft = log(abs(fftshift(fftn(image_fft))));

    [Height, Width, Depth] = size(image_fft);
    minval = min(image_fft, [], 'all');
    maxval = max(image_fft, [], 'all');
    image_fft = permute(image_fft, [2 1 3]);
    image_fft = image_fft(:);
    image_fft = reshape(image_fft, [1, size(image_fft)]);
    
    fileID = fopen(strcat('Wiener_dir_', num2str(dir), '_fft.tif'), 'w+');
    header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(image_fft), ...
        1, Depth, 1, imgParams.dz / params.z_zoom, minval, maxval, imgParams.dy / params.zoomfact);
    
    [~, ~, system_endian] = computer;
    if system_endian ~= header.endian
        image_fft = swapbytes(image_fft);
    end
    
    image_fft = typecast(image_fft, 'uint8');
    fwrite(fileID, image_fft, 'uint8');
    fclose(fileID);
end

% Wrtie per-direction Wiener filtered image
if ~(params.bkeepnegative)
    image(image < 0) = 0;
end

[Height, Width, Depth] = size(image);
minval = min(image, [], 'all');
maxval = max(image, [], 'all');
image = permute(image, [2 1 3]);
image = image(:);
image = reshape(image, [1, size(image)]);

fileID = fopen(strcat(params.ifiles, '\Wiener_dir_', num2str(dir), '_', num2str(it), '.tif'), 'w+');
header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(image), ...
    1, Depth, 1, imgParams.dz / params.z_zoom, minval, maxval, imgParams.dy / params.zoomfact);

[~, ~, system_endian] = computer;
if system_endian ~= header.endian
    image = swapbytes(image);
end

image = typecast(image, 'uint8');
fwrite(fileID, image, 'uint8');
fclose(fileID);

end

%% Save per-direction Wiener filtered result
function saveWienerResult_single_direction_shifted(dir, it, infft, params, imgParams)

image = real(ifftn(infft));

% Write per-direction Wiener filtered fft data image
if (params.bsavefft)
    if ~(params.bfftkeepnegative)
        infft = image;
        infft(infft < 0) = 0;
        infft = log(abs(fftshift(fftn(infft))));
    else
        infft = log(abs(fftshift(infft)));
    end
    
    [Height, Width, Depth] = size(infft);
    minval = min(infft, [], 'all');
    maxval = max(infft, [], 'all');
    infft = permute(infft, [2 1 3]);
    infft = infft(:);
    infft = reshape(infft, [1, size(infft)]);
    
    fileID = fopen(strcat('Wiener_dir_', num2str(dir), '_fft.tif'), 'w+');
    header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(infft), ...
        1, Depth, 1, imgParams.dz / params.z_zoom, minval, maxval, imgParams.dy / params.zoomfact);
    
    [~, ~, system_endian] = computer;
    if system_endian ~= header.endian
        infft = swapbytes(infft);
    end
    
    infft = typecast(infft, 'uint8');
    fwrite(fileID, infft, 'uint8');
    fclose(fileID);
end

% Write per-direction Wiener filtered image
if ~(params.bkeepnegative)
    image(image < 0) = 0;
end

[Height, Width, Depth] = size(image);
minval = min(image, [], 'all');
maxval = max(image, [], 'all');
image = permute(image, [2 1 3]);
image = image(:);
image = reshape(image, [1, size(image)]);

fileID = fopen(strcat(params.ifiles, '\Wiener_dir_', num2str(dir), '_', num2str(it), '.tif'), 'w+');
header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(image), ...
    1, Depth, 1, imgParams.dz / params.z_zoom, minval, maxval, imgParams.dy / params.zoomfact);

[~, ~, system_endian] = computer;
if system_endian ~= header.endian
    image = swapbytes(image);
end

image = typecast(image, 'uint8');
fwrite(fileID, image, 'uint8');
fclose(fileID);

end

%% to pre-filter separated bands for all directions;
function bands = filterbands(dir, bands, k0, amp, noiseVarFactors, otf, params, imgParams) 
%{
    dir = 1, 2, 3; Suppose dir = 1;
    bands = obj.m_reconData.savedBands(:, :, :, :, direction);
    k0 =  obj.m_reconData.k0;
    amp = obj.m_reconData.amp;
    noiseVarFactors = obj.m_reconData.noiseVarFactors;
    otf = obj.otf;
    params = obj.m_myParams;
    imgParams = obj.m_imgParams;
%}

ampmag2 = zeros([params.norders, 1], 'single');
conjamp = complex(zeros([params.norders, 1], 'single'), zeros([params.norders, 1], 'single'));

for order = 1 : params.norders
    ampmag2(order) = abs(amp(order, dir)) ^ 2;
    conjamp(order) = conj(amp(order, dir));
end

% Explicitly calculate mag2 of amp for all orders
ampmag2_alldirs = abs(amp) .^ 2;

for order = 0 : (params.norders - 1)
    fprintf("dir %d: Wiener filtering order %d\n", dir, order);
    
    if (order == 0)
        dev_bandptr = bands(:, :, :, 1);
        dev_bandptr2 = complex(zeros(size(dev_bandptr), 'single'), zeros(size(dev_bandptr), 'single'));
        
        [dev_bandptr, ~] = filterbands_kernel1(order, dir, k0, dev_bandptr, dev_bandptr2, ...
            ampmag2, noiseVarFactors, ampmag2_alldirs, conjamp, otf, params, imgParams);
        
        if (imgParams.nz0 > (2 * imgParams.zdistcutoff(order + 1) + 1))
            [dev_bandptr, ~] = filterbands_kernel3(order, dev_bandptr, dev_bandptr2, imgParams.zdistcutoff, imgParams.nz0);
        end
        
        bands(:, :, :, 1) = dev_bandptr;
    else
        dev_bandptr = bands(:, :, :, 2 * order);
        dev_bandptr2 = bands(:, :, :, 2 * order + 1);
        
        [dev_bandptr, dev_bandptr2] = filterbands_kernel1(order, dir, k0, dev_bandptr, dev_bandptr2, ...
            ampmag2, noiseVarFactors, ampmag2_alldirs, conjamp, otf, params, imgParams);
        
        if (imgParams.nz0 > (2 * imgParams.zdistcutoff(order + 1) + 1))
            [dev_bandptr, dev_bandptr2] = filterbands_kernel3(order, dev_bandptr, dev_bandptr2, imgParams.zdistcutoff, imgParams.nz0);
        end
        
        bands(:, :, :, 2 * order) = dev_bandptr;
        bands(:, :, :, 2 * order + 1) = dev_bandptr2;
    end
end

end

%% to compute wiener filter components
function [band, band2] = filterbands_kernel1(order, dir, k0, band, band2, ampmag2, noiseVarFactors, ampmag2_alldirs, conjamp, otf, params, imgParams)
%{
    k0 = obj.m_reconData.k0;
    dir = 1, 2, 3; Suppose dir = 1;
    ndirs = obj.m_myParams.ndirs = 3;
    order = 0, 1, 2;

    bSecondEntry = false (1st call) or true (2nd call);
    OTF = obj.m_reconData.otf(:, :, :, :, dir_);
    ampmag2: magitude square of amp for all orders
    noiseVarFactors = obj.m_reconData.noiseVarFactors;
    params = obj.m_myParams;
    imgParams = obj.m_imgParams;
%}

dir_otf = 1;
if (params.bOneOTFperAngle)
    dir_otf = dir;
end

kx = order * k0(dir).x;
ky = order * k0(dir).y;

% otf1 is the order = 0, 1, 2 OTF value corresponding to the current (x1, y1, z0) or say "relative index to center point (0, 0, 0)" in image;
otf1 = otf.getOtfMatrix3D(order, dir_otf, imgParams.nx, imgParams.ny, 0, 0, imgParams.krscale, imgParams.zdistcutoff(order + 1), imgParams.kzscale);
% otf1: OTF(order = m)(k)
weight = abs(otf1) .^ 2;
if (order ~= 0)
    weight = weight .* ampmag2(order + 1);
end

dampen_matrix = otf.getAttenuationMatrix3D(order, imgParams.nx, imgParams.ny, 0, 0, imgParams.zdistcutoff(order + 1), imgParams.rdistcutoff, params);
dampfact = dampen_matrix / noiseVarFactors(dir, order + 1);

weight = weight .* dampfact;
sumweight = weight;

% Calculate all-direction Wiener denominator
for dir2 = 1 : params.ndirs
    for order2 = -(params.norders - 1) : (params.norders - 1)
        % order2 = [-2, -1, 0, 1, 2];
        if (dir2 == dir && order2 == order) 
            continue;
        end
        
        %  bFilteroverlaps is always true except when (during debug) generating an unfiltered exploded view 
        if (~params.bFilteroverlaps && ~(order2 == 0 && order == 0))
            continue;
        end
        
        % shifted center of band 2
        kx2 = order2 * k0(dir2).x;
        ky2 = order2 * k0(dir2).y;
        
        dir_otf = 1;
        if (params.bOneOTFperAngle)
            dir_otf = dir2;
        end
        
        % otf2 is the abs(order2) = 0, 1, 2 OTF value corresponding to the current (x2, y2, z0) or say "relative index to center point (0, 0, 0)" in image;
        otf2 = otf.getOtfMatrix3D(abs(order2), dir_otf, imgParams.nx, imgParams.ny, kx2 - kx, ky2 - ky, imgParams.krscale, imgParams.zdistcutoff(order + 1), imgParams.kzscale);
        weight = abs(otf2) .^ 2;
        weight = weight ./ noiseVarFactors(dir2, abs(order2) + 1);
        
        amp2mag2 = ampmag2_alldirs(abs(order2) + 1, dir2);
        if (order2 ~= 0)
            weight = weight .* amp2mag2;
        end
        
        dampen_matrix = otf.getAttenuationMatrix3D(abs(order2), imgParams.nx, imgParams.ny, kx2 - kx, ky2 - ky, imgParams.zdistcutoff(order + 1), imgParams.rdistcutoff, params);
        weight = weight .* dampen_matrix;
        sumweight = sumweight + weight;
    end
end
% After the two for loops: sumweight = dampfact * (otf1.x ^ 2 + otf1.y ^ 2) * Magnitude[m_reconData.amp[dir][order]] + Sum[dev_suppress(rdist2) * (otf2.x ^ 2 + otf2.y ^ 2) / 1.0f * Magnitude[m_reconData.amp[dir2 = 0,1,2][order2 = 0,1,2]]];
% 3 directions * 5 orders = 15 terms in total;

sumweight = sumweight + params.wiener ^ 2;
% Above is the denominator part of the shifted wiener filter function;

scale = dampfact .* conj(otf1);
% scale = dampfact .* conj(OTF(dir, order) = dampfact .* O(*)(m)(k)
scale = scale ./ sumweight;
% Based on 2008 Biophysical Journal paper's defination:
% scale = conj[OTF(dir, order)] / (sum[abs(OTF(dir2, order2) shifted by [-order2 * p(dir2) + order * k0(dir)]) .^ 2] + wiener ^ 2);


% Apodization function: A(k - mp(d)) = A(k - order * k0), which is function shifted by -mp(d)
if (params.apodizeoutput)
    apofact = otf.writeApoFactor(imgParams.nx, imgParams.ny, -kx, -ky, imgParams.zdistcutoff(order + 1), imgParams.apocutoff, imgParams.zapocutoff, params);
    scale = scale .* apofact;
    % Based on 2008 Biophysical Journal paper's defination:
    % scale = conj[OTF(order)] * A(k - order * p(dir)) /
    % (sum[abs(OTF(dir2, order2) shifted by [order2 * p(dir2) - order * k0(dir)]) .^ 2] + wiener ^ 2);
    % Simply speaking, scale = A(k - order * p(dir)) * conjugate(otf1)/sumweight
end
   
z0 = -imgParams.zdistcutoff(order + 1) : imgParams.zdistcutoff(order + 1);  % z0 = -20:0:20
iz = mod((z0 + imgParams.nz0), imgParams.nz0);          % z0 = [-20:20] --> [80:99]...[0:20] = iz

%{
    separate (for this pixel) the even and odd "bands" into the true plus and minus bands
	apply the scale to the plus band only (this will apply it to the minus band as well by symmetry(?)
    reassemble into the even and odd "bands"
%}
if order == 0
    band(:, :, iz + 1) = band(:, :, iz + 1) .* scale;
else
    scale = scale .* conjamp(order + 1);
    % Based on 2008 Biophysical Journal paper's defination:
    % scale = A(k - order * p(dir)) * [c(m)* exp(-1j * phi(m))] * conjugate(otf1)/sumweight
    % conjamp(order + 1) = c(m)* exp(1j * phi(m))
    
    % Get the "genuine" seprated plus bands
    bandplus = band - 1j .* band2;
    
    % scale only the bandplus part
    bandplus(:, :, iz + 1) = bandplus(:, :, iz + 1) .* scale;
    band2 = bandplus;
end

end

%% to clear everything above and below zdistcutoff to 0
function [band, band2] = filterbands_kernel3(order, band, band2, zdistcutoff, nz)

z0 = (1 + zdistcutoff(order + 1)) : (nz - zdistcutoff(order + 1) -  1) ;
band(:, :, z0 + 1) = 0;
if (order ~= 0)
    band2(:, :, z0 + 1) = 0;
end

end

%% to assemble the bands together;
function outbuffer = assemblerealspacebands(dir, outbuffer, bigbuffer, bands, k0, params, imgParams)
%{
    dir = 1, 2, 3;
	outbuffer = obj.m_reconData.outbuffer = zeros([2 * 512, 2 * 512, 1 * 100], 'single');
	bigbuffer = obj.m_reconData.bigbuffer = complex(zeros([2 * 512, 2 * 512, 1 * 100], 'single'), ...);
	bands = obj.m_reconData.savedBands(:, :, :, :, direction);
	k0 = obj.m_reconData.k0;

    params = obj.m_myParams;
    imgParams = obj.m_imgParams;
%}

    % expfact is used for "exploded view".  For normal reconstruction expfact = 1.0
	fact = params.explodefact / 0.5;
    
    % 0th order
    fprintf("moving centerband\n");
    bigbuffer(:) = 0;
    bigbuffer = placeFreq(bands(:, :, :, 1), bigbuffer);

    % transform it into real space
	fprintf("re-transforming central band\n");
    bigbuffer = ifftn(bigbuffer);
    
    fprintf("inserting central band\n");
    outbuffer = write_outbuffer_kernel1(bigbuffer, outbuffer);
    
    fprintf("central band assembly completed\n");
    
    % Side bands
    for order = 1 : (params.norders - 1)
        fprintf("moving order %d\n", order);
        bigbuffer(:) = 0;
        bigbuffer = placeFreq(bands(:, :, :, 2 * order + 1), bigbuffer);
        
        % transform it into real space
        fprintf("re-transforming side bands\n");
        bigbuffer = ifftn(bigbuffer);
        
        % For 3D, prepare 2D array of sines and cosines first, then loop over z.
        k0x = k0(dir).x * order;
        k0y = k0(dir).y * order;
        x = 0 : (params.zoomfact * imgParams.nx - 1);
        y = 0 : (params.zoomfact * imgParams.ny - 1);
        z = 0 : (params.z_zoom * imgParams.nz0 - 1);
        [X, Y, ~] = meshgrid(x, y, z);
        angle = fact * pi * (X .* k0x / (params.zoomfact * imgParams.nx) + Y .* k0y / (params.zoomfact * imgParams.ny));
        coslookup = cos(-angle);
        sinlookup = sin(-angle);
        
        fprintf("inserting side bands\n");
        outbuffer = write_outbuffer_kernel2(coslookup, sinlookup, bigbuffer, outbuffer);
        
        fprintf("order %d side band assembly completed\n", order);
    end
    
    fprintf("\n");

end

%% to add the real part of m_reconData.bigbuffer (0th order) into m_reconData.outbuffer; */
function outbuffer = write_outbuffer_kernel1(bigbuffer, outbuffer)

outbuffer = outbuffer + real(bigbuffer);

end

%% to add the real part of m_reconData.bigbuffer (1st & 2nd ... orders) into m_reconData.outbuffer; */
function outbuffer  = write_outbuffer_kernel2(coslookup, sinlookup, bigbuffer, outbuffer)

outbuffer = outbuffer + 2 * real(bigbuffer .* (coslookup + 1j * sinlookup));

end

%% to use 'fullResult' as numerator and calculate (1) Wiener denominator & (2) Apodization matrix to do Wiener filter
function outbuffer = filterbands_shifted(outbuffer, fullResult, k0, amp, noiseVarFactors, otf, params, imgParams)

% Calculate all-direction Wiener denominator
denom = params.wiener ^ 2;

for dir = 1 : params.ndirs
    dir_otf = 1;
    if (params.bOneOTFperAngle)
        dir_otf = dir;
    end
    
    for order = -(params.norders - 1) : (params.norders - 1)
        kx = order * k0(dir).x;
        ky = order * k0(dir).y;
        
        OTF = zeros(size(outbuffer), 'single');
        otfval = otf.getOtfMatrix3D(abs(order), dir_otf, params.zoomfact * imgParams.nx, params.zoomfact * imgParams.ny, kx, ky, imgParams.krscale, imgParams.zdistcutoff(abs(order) + 1), imgParams.kzscale);
        otfval = abs(otfval) .^2 ./ noiseVarFactors(dir, abs(order) + 1);
        if (order ~= 0)
            otfval = otfval .* abs(amp(abs(order) + 1, dir)) ^ 2;
        end
        
        dampen_matrix = otf.getAttenuationMatrix3D(abs(order), params.zoomfact * imgParams.nx, params.zoomfact * imgParams.ny, kx, ky, imgParams.zdistcutoff(abs(order) + 1), imgParams.rdistcutoff, params);
        otfval = otfval .* dampen_matrix;
        
        z0 = -imgParams.zdistcutoff(abs(order) + 1) : imgParams.zdistcutoff(abs(order) + 1);
        iz = mod((z0 + params.z_zoom * imgParams.nz0), params.z_zoom * imgParams.nz0);
        OTF(:, :, iz + 1) = otfval;
        denom = denom + OTF;
    end
end

outbuffer = fullResult ./ denom;

if (params.apodizeoutput)
    apo = zeros(size(outbuffer), 'single');
    apofact = otf.writeApoFactor(params.zoomfact * imgParams.nx, params.zoomfact * imgParams.ny, 0, 0, imgParams.zapocutoff, imgParams.apocutoff, imgParams.zapocutoff, params);
    
    z0 = -imgParams.zapocutoff : imgParams.zapocutoff;
    iz = mod((z0 + params.z_zoom * imgParams.nz0), params.z_zoom * imgParams.nz0);
    apo(:, :, iz + 1) = apofact;
    
    outbuffer = outbuffer .* apo;
end

end

%% Run Richardson-Lucy deconvolution iter_num on img. 
function estimate = deconvolve(estimate, otf, iter_num, inputIsInFreqSpace)
%{
    estimate:   Input image, will be modified and contain result, could be 2D or 3D
    otf:        The optical transfer function to use, must have the same dimension as img
    iter_num:   Maximum number of iteration steps
    inputIsInFreqSpace: If input is already in freq. space, will set output to same space
%}

if (size(estimate) ~= size(otf))
    error('In RL deconvolution, the size of input image in not equal to otf.');
end

Ratio = ones(size(estimate), 'single');
Normalization = abs(ifftn(fftn(Ratio) .* conj(otf)));

if (inputIsInFreqSpace)
    estimate = abs(ifftn(estimate));
end

for iteration = 1 : iter_num
    if (iter_num > 1)
        fprintf('   processing iteration # %d...\n', iteration);
    end
    
    % 1st step: compute Blur = Estimate (*) PSF
    Blur = abs(ifftn(fftn(estimate) .* otf));
    
    % 2nd step: compute Ratio = Image / Blur
    Ratio = estimate ./ Blur;
    
    % 3rd step: compute Correction = Ratio (*) PSF_flip
    Correction = abs(ifftn(fftn(Ratio) .* conj(otf)));
    
    % 4th step: compute Next_Estimate = Estimate .* Correction ./ Normalization
    estimate = estimate .* Correction ./ (Normalization + realmin('single'));
    
    estimate(estimate < 0) = realmin('single');
end

end
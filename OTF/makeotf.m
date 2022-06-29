% Matlab version of Mats and Lin's 'makeotf.c'

% A new small project to radially average OTF starting from PSF dataset.
% Mainly for structured illumination microscopy (1 or 2 objectives).
% Also take into account the following:
%   1.  compensation for finite bead size;
%   2.  the phase factor associated with the side bands;
%   3.  estimate the sub-pixel position of the bead by fitting parabolas;
%   4.  cleanup out-of-band noises.

% clear all;
close all;
clc;

% NA = 1.35;          % Emission numerical aperture
NA = 1.27;
lambda = 488;        % Excitation wavelength
% lamda = 561;
twolens = false;    % Whether using mirror-based 4-beam SIM acquisition

%% Pre-defined constants
% cd D:\Matlab_Code\1P27NA_water\OTF;
MAXDIRs = 3;
SPOTRATIO = 3;
MAXPHASES = 10;
DEFAULTPHASES = 5;

%% Variables
ifiles = []; ofiles = []; I2Mfiles = []; bandsfiles = [];
istream_no = []; ostream_no = [];
nx = []; ny = []; nz = []; nxy = []; phase = []; nphases = DEFAULTPHASES; norders = []; order = []; napodize = 10;
i = []; j = []; z = []; border_size = [];
dr = []; dz = []; dkr = []; dkz = []; background = []; background_dummy = -1.0; background_i2m = [];
xcofm = []; ycofm = []; zcofm = []; xcofm_i2m = []; ycofm_i2m = []; zcofm_i2m = [];
floatimage = []; sepMatrix = []; floatsection = []; I2M_image = [];
buffer = [];
bands = []; avg_output = []; I2Mavg_output = [];
bead_diameter = 0.100;                  % Bead diameter is 0.100 um (100 nm);
scalefactor = 1.0;

if NA == 1.27
    k0angleguess = (-80.2723) / 180 * pi;   % Estimated angle of the 1st orientation
    NIMM = 1.333;                       % Refractive index of immersion medium
    if lambda == 488
        wavelength = 525;               % Fluorescence emmision wavelength 488 nm / 525 nm;
        linespacing = 0.209169;         % Line spacing of SIM pattern in microns, 0.20911 um for 488 nm channel 1.27 NA.
        if twolens
            ifiles = 'PSF_4BSIM_1P27NA_488.tif';
            ofiles = 'OTF_4BSIM_1P27NA_488.tif';
            ileavekz = [43, 44, 11, 12, 32, 33, 0]; % wavelength = 525; Water (RI = 1.333);
        else
            ifiles = 'PSF_3DSIM_1P27NA_488.tif';
            ofiles = 'OTF_3DSIM_1P27NA_488.tif';
            ileavekz = [0, 0, 11, 12, 0, 0, 0];     % wavelength = 525; Water (RI = 1.333)
        end
    elseif lambda == 561
        wavelength = 607;               % Fluorescence excitation/emmision wavelength 561 nm / 607 nm;
        linespacing = 0.240117;         % Line spacing of SIM pattern in microns, 0.240117 um for 561 nm channel 1.27 NA.
        if twolens
            ifiles = 'PSF_4BSIM_1P27NA_561.tif';
            ofiles = 'OTF_4BSIM_1P27NA_561.tif';
            ileavekz = [37, 38, 9, 10, 28, 29, 0];  % wavelength = 607; Water (RI = 1.333);
        else
            ifiles = 'PSF_3DSIM_1P27NA_561.tif';
            ofiles = 'OTF_3DSIM_1P27NA_561.tif';
            ileavekz = [0, 0, 9, 10, 0, 0, 0];      % wavelength = 607; Water (RI = 1.333)
        end
    else
        error('No PSF data with such excitation wavelength');
    end
elseif NA == 1.35
    k0angleguess = (-11.062) / 180 * pi;    % Estimated angle of the 1st orientation
    NIMM = 1.406;                       % Refractive index of immersion medium
    if lambda == 488
        wavelength = 525;               % Fluorescence emmision wavelength 488 nm / 525 nm;
        linespacing = 0.199378;         % Line spacing of SIM pattern in microns, 0.199315 for 488 nm channel 1.35 NA.
        if twolens
            ifiles = 'PSF_4BSIM_1P35NA_488.tif';
            ofiles = 'OTF_4BSIM_1P35NA_488.tif';
            ileavekz = [47, 48, 12, 13, 35, 36, 0]; % wavelength = 525; Iodixanol (45.6% w/v RI = 1.406)
        else
            ifiles = 'PSF_3DSIM_1P35NA_488.tif';
            ofiles = 'OTF_3DSIM_1P35NA_488.tif';
            ileavekz = [0, 0, 12, 13, 0, 0, 0];     % wavelength = 525; Iodixanol (45.6% w/v RI = 1.406)
        end
    elseif lambda == 561
        wavelength = 607;               % Fluorescence excitation/emmision wavelength 561 nm / 607 nm;
        linespacing = 0.228782;         % Line spacing of SIM pattern in microns, 0.228782 for 561 nm channel 1.35 NA.
        if twolens
            ifiles = 'PSF_4BSIM_1P35NA_561.tif';
            ofiles = 'OTF_4BSIM_1P35NA_561.tif';
            ileavekz = [40, 42, 10, 11, 30, 31, 0]; % wavelength = 607; Iodixanol (45.6% w/v RI = 1.406)
        else
            ifiles = 'PSF_3DSIM_1P35NA_561.tif';
            ofiles = 'OTF_3DSIM_1P35NA_561.tif';
            ileavekz = [0, 0, 10, 11, 0, 0, 0];     % wavelength = 607; Iodixanol (45.6% w/v RI = 1.406)
        end
    else
        error('No PSF data with such excitation wavelength');
    end
else
    error('No PSF data with such NA');
end

% ileavekz = [0, 0, 0, 0, 0, 0, 0];


header = Simple_IFD(); otfheader = [];

ifixz = []; ifixr = [];
dorescale = 0;                          % Whether rescale the OTF for order0, order1 and order2 respectively
five_bands = 0;                         % Whether to create 5-phase bands or 3-order bands OTF
do_compen = 1;                          % Whether to compensate finite bead size
I2M_inc = 0;
Generate_bands = true;                  % Whether to save the band-separated input PSF image
conjugate = 0;                          % Wheter to calculate the conjugatation of OTF
bBgInExtHdr = 0;                        % If the background of each section is recorded in the extended header's 3rd float number (in Lin's scope)
zsec = [];

icleanup = [];
ileave1_1 = []; ileave1_2 = []; ileave2 = [];
%{
    ileave1_1 = ileavekz[0] is the minimum shifted-pixel number for order1 OTF;
    ileave1_2 = ileavekz[1] is the maximum shifted-pixel number for order1 OTF;
    ileave2 = ileavekz[2] is the maximum axially-broadened pixel number for order2 OTF;
%}
interpkr = [];                          % Whether to fix origin

bUseCorr = 0;                           % Whether to use a camera flat-fielding (or correction) file, default = 0
bForcedPIshift = 0;
corrfiles = [];
background2D = 0;
slope2D = 0;

interpkr = [0, 0];

%% Variable assignment
interpkr(1) = 3;
interpkr(2) = 20;
ifixz(1) = 	1; ifixz(2) = 1; ifixz(3) = 1; ifixr = 0;

ofiles_exp = 'OTF_exp.tif';
bandsfiles = 'bands.tif';


%% PSF and OTF dataset file reading / creating.
if (nphases > MAXPHASES)
    warning('nphases is larger than MAXPHASES');
    return;
end

if not(isempty(ifiles))
    fprintf('PSF dataset file name: %s\n', ifiles);
end

[istream_no, errmsg] = fopen(ifiles, 'r');
if not(isempty(errmsg))
    warning(errmsg);
    warning('File: %s does not exist.', ifiles);
    return;
end

if not(isempty(ofiles))
    fprintf('Output OTF file name: %s\n', ofiles);
end

[ostream_no, errmsg] = fopen(ofiles, 'w+');
if not(isempty(errmsg))
    warning(errmsg);
    warning('File: %s can not be created.', ofiles);
    return;
end

if not(isempty(ofiles_exp))
    fprintf('Output OTF exponential file name: %s\n', ofiles);
end

[ostream_exp_no, errmsg] = fopen(ofiles_exp, 'w+');
if not(isempty(errmsg))
    warning(errmsg);
    warning('File: %s can not be created.', ofiles_exp);
    return;
end

%% Get the entire header in a single block of memory
header = ImageJ_formatted_TIFF.parse_tif(istream_no, 0);

nx = header.ImageWidth;
ny = header.ImageLength;
ImageDescription = convertStringsToChars(header.ImageDescription);
% Calculate nz from header.ImageDescription
k1 = strfind(ImageDescription, 'images=');
ImageDescription_crop = ImageDescription(k1:end);
k2 = strfind(ImageDescription_crop, newline);
if ~isempty(k1) && ~isempty(k2)
    nz = str2double(ImageDescription(k1(1) + 7: k1(1) + k2(1) - 2));
else
    warning("Did not find 'images=' in ImageDescription. Try to calculate it from filesize.");
    FileAttributes = dir(ifiles);
    nz = round(FileAttributes.bytes / (header.ImageWidth * header.ImageLength * (header.BitsPerSample / 8)));
end
% Calculate spacing from header.
k1 = strfind(ImageDescription, 'spacing=');
ImageDescription_crop = ImageDescription(k1:end);
k2 = strfind(ImageDescription_crop, newline);
if ~isempty(k1) && ~isempty(k2)
    spacing = single(str2double(ImageDescription(k1(1) + 8: k1(1) + k2(1) - 2)));
else
    warning("Did not find 'spacing=' in ImageDescription. Suppose it is 0.04 microns.");
    spacing = single(0.04);
end
% Calculate data_type from header.BitsPerSample
data_type = [header.SampleFormat, header.BitsPerSample];
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
    error("ReadTifStack does not support SampleFormat:%d with BitsPerSample: %d", header.SampleFormat, header.BitsPerSample);
end

% Initialize reading buffer and parameters
PixelNum = header.ImageWidth * header.ImageLength;
ByteCounts = fix(PixelNum * header.BitsPerSample / 8);
fseek(istream_no, header.StripOffsets, 'bof');
[~, ~, system_endian] = computer;


if ~I2M_inc
    nz = nz / nphases;              % nz = 320 / 5;
else
    nz = nz / (nphases + 1);        % nz = 384 / 5;
end

dr = header.resolution;
dz = spacing;

dkr = 1 / (ny * dr);                % dkr = 1 / (256 * 0.08) = 1 / 20.48;
dkz = 1 / (nz * dz);                % dkz = 1 / (64 * 0.125) = 1 / 8;

fprintf("nx=%d, ny=%d, nz=%d\n", nx, ny, nz);

norders = (nphases + 1) / 2;        % norders = (5 + 1) / 2 = 3;
sepMatrix = makematrix(nphases);
nxy = (nx + 2) * ny;                % nxy = (256 + 2) * 256;



%% Allocate memory for floatimage and bands
buffer = zeros([ny, nx], 'single');
floatimage = zeros([ny, nx, nz, nphases], 'single');
bands = complex(zeros(size(floatimage), 'single'));
if I2M_inc
    I2M_image = zeros([ny, nx, nz], 'single');
end

%% Allocate memory for background
background = zeros([nz, nphases], 'single');
border_size = fix(20 / 256 * ny);   % border_size = fix(20 / 256 * 256) = 20;

%% Flatfield correction of measured data using calibration data
if bUseCorr
    disp('Loading sCMOS calibration file.\n');
    [background2D, slope2D] = getbg_and_slope(corrfiles, nx, ny);
end

fprintf("Reading data...\n")
napodize = fix(napodize / 256 * ny);
zsec = 0;
for z = 1:nz
    for phase = 1:nphases
        % phase first, then z; data organized into (nz, ndirs, nphases)
        
        % Reads the next section into ImgBuffer and advances the file pointer to the section after that.
        buffer_uint16 = typecast(transpose(fread(istream_no, ByteCounts, 'uint8=>uint8')), dtype);
        if system_endian ~= header.endian
            buffer_uint16 = swapbytes(buffer_uint16);
        end
        buffer(:) = single(buffer_uint16);
        
        if bBgInExtHdr
        elseif bUseCorr
            buffer = buffer - background2D;
            buffer(buffer < 0) = 0;
            buffer = buffer .* slope2D;     % buffer = (buffer - background2D) * slope2D;
        elseif background_dummy >= 0
            background(z, phase) = background_dummy;    % The artifical background photon numbers we assumed;
        else
            % background is estimated using the average intensities of the 20 marginal pixels
            background(z, phase) = estimate_background(buffer, nx, ny, border_size);
        end
        
        if napodize >= 0
            floatsection = apodize(napodize, nx, ny, buffer);
        elseif napodize == -1
            floatsection = cosapodize(nx, ny, buffer);
        end
        
        floatimage(:, :, z, phase) = floatsection;
        zsec = zsec + 1;
    end
    
    % Data contains one extra section of I2M image
    if I2M_inc
        buffer_uint16 = typecast(transpose(fread(istream_no, ByteCounts, 'uint8=>uint8')), dtype);
        if system_endian ~= header.endian
            buffer_uint16 = swapbytes(buffer_uint16);
        end
        buffer(:) = single(buffer_uint16);
        I2M_image(:, :, z) = buffer;
    end
end
fclose(istream_no);

%% Before FFT, use center band to estimate bead center position
% you can choose either "parabolic fitting" or "window binning maximum" method, only very little difference
% [xcofm, ycofm, zcofm, background_dummy, background_i2m, xcofm_i2m, ycofm_i2m, zcofm_i2m] = determine_center_and_background(floatimage, I2M_image, nx, ny, nz, nphases, twolens, I2M_inc);
[xcofm, ycofm, zcofm, background_dummy, background_i2m, xcofm_i2m, ycofm_i2m, zcofm_i2m] = determine_center_and_background_2(floatimage, I2M_image, nx, ny, nz, nphases, I2M_inc);

fprintf("Center of mass is (%.3f, %.3f, %.3f)\n", xcofm, ycofm, zcofm);
fprintf("Background is %.3f\n", background_dummy);

if I2M_inc
    fprintf("I2M psf's background is %.3f\n", background_i2m);
    fprintf("I2M psf's center of mass is (%.3f, %.3f, %.3f)\n\n", xcofm_i2m, ycofm_i2m, zcofm_i2m);
end

%% Subtract background and seperate bands
for z = 1:nz
    for phase = 1:nphases
        floatimage(:, :, z, phase) = floatimage(:, :, z, phase) - background(z, phase);
    end
    
    if I2M_inc
        I2M_image = I2M_image - backgournd_i2m;
        I2M_image(I2M_image < 0) = 0;
    end
end

floatimage(floatimage < 0) = 0;
if I2M_inc
    I2M_image(I2M_image < 0) = 0;
end

%% Mathematically, band separation in real vs. Fourier space is no difference
% bandSepInFFtSpace = false;
% if (bandSepInFFtSpace)
%     for phase = 1 : nphases
%         bands(:, :, :, phase) = fftn(floatimage(:, :, :, phase));
%     end
% end

for z = 1 : nz
    % After this function, we separate bands at each z layer
    floatimage(:, :, z, :) = separate(nx, ny, z, nphases, floatimage, sepMatrix);
end

% if (bandSepInFFtSpace)
%     for phase = 1 : nphases
%         floatimage(:, :, :, phase) = ifftn(bands(:, :, :, phase));
%     end
% end


%% Just write bands 0, 1 & 2 into file
if (Generate_bands)
    fileID = fopen(bandsfiles, 'w+');
    minval = min(abs(floatimage), [], 'all');
    maxval = max(abs(floatimage), [], 'all');
    header_bands = ImageJ_formatted_TIFF.write_IFD(fileID, nx, ny, class(floatimage), 3, nz, 1, dz, minval, maxval, dr);
    Stack_bands = permute(floatimage(:, :, :, [1 2 4]), [2 1 4 3]);
    Stack_bands = reshape(Stack_bands, [1, nx * ny * nz * 3]);
    
    if system_endian ~= header_bands.endian
        Stack_bands = swapbytes(Stack_bands);
    end
    Stack_bands = typecast(Stack_bands, 'uint8');
    fwrite(fileID, Stack_bands, 'uint8');
    fclose(fileID);
end

%% 3D-FFT from floatimage into bands
for phase = 1 : nphases
    bands(:, :, :, phase) = fftn(floatimage(:, :, :, phase));
end

if I2M_inc
    I2M_image = fftn(I2M_image);
end

%% Modify the phase of bands, so that it corresponds to FFT of a bead at origin
fprintf("Shifting center...\n");
for phase = 1 : nphases
    bands(:, :, :, phase) = shift_center(bands(:, :, :, phase), nx, ny, nz, xcofm, ycofm, zcofm);
end

if I2M_inc
    I2M_image = shift_center(I2M_image, nx, ny, nz, xcofm_i2m, ycofm_i2m, zcofm_i2m);
end

%% Now compensate finite bead size for all bands
if do_compen
    bands = beadsize_compensate(bands, k0angleguess, linespacing, bead_diameter, norders, nx, ny, nz, dkr, dkz);
end

%% Radial average
avg_output = complex(zeros([fix(nx / 2) + 1, nz, nphases], 'single'));  % [129, 64, 5]
for order = 0: (norders - 1)
    if order == 0
        avg_output(:, :, 1) = radialft(bands(:, :, :, 1), nx, ny, nz);
    else
        avg_output(:, :, 2 * order) = radialft(bands(:, :, :, 2 * order), nx, ny, nz);
        avg_output(:, :, 2 * order + 1) = radialft(bands(:, :, :, 2 * order + 1), nx, ny, nz);
    end
end

if I2M_inc
    I2Mavg_output = complex(zeros([fix(nx / 2) + 1, nz], 'single'));   % [129, 64, 5]
    I2Mavg_output(:) = radialft(I2M_image, nx, ny, nz);
end


%% Modify, cleanup, fixorigin (optional) and rescale ;
icleanup = fix(nx / 2) + 1;     % icleanup = 129;
%{
    ileave1_1 = ileavekz(1) is the minimum shifted-pixel number for order1 OTF;
    ileave1_2 = ileavekz(2) is the maximum shifted-pixel number for order1 OTF;
    ileave2 = ileavekz(3) is the maximum axially-broadened pixel number for order2 OTF;
%}

ileave0_1 = ileavekz(1);
ileave0_2 = ileavekz(2);
ileave1_1 = ileavekz(3);
ileave1_2 = ileavekz(4);
ileave1_3 = ileavekz(5);
ileave1_4 = ileavekz(6);
ileave2 = ileavekz(7);

for phase = 1 : nphases
    avg_output(:, :, phase) = modify(avg_output(:, :, phase), nx, nz, ifixz, ifixr, fix(phase / 2), twolens);
    
    if (ileave1_1 > 0)
        avg_output(:, :, phase) = cleanup(avg_output(:, :, phase), fix(phase / 2), nx, nz, dkr, dkz,...
            linespacing, wavelength, icleanup, ileave0_1, ileave0_2, ileave1_1, ileave1_2, ileave1_3, ileave1_4, ileave2, twolens, NA, NIMM);
    end
    
    if (phase == 1) && (interpkr(1) > 0)
        avg_output(:, :, phase) = fixorigin(avg_output(:, :, phase), nx, nz, interpkr(1), interpkr(2), dkr, dkz, wavelength, twolens, NA, NIMM);
    end
    
    [avg_output(:, :, phase), scalefactor] = rescale(avg_output(:, :, phase), fix(phase / 2), nx, nz, scalefactor, dorescale);
end

%% For side bands, combine bandre's (real band) and bandim's (imaginary band) into bandplus
if (~five_bands)
    avg_output = combine_reim(avg_output, norders, nx, nz, bForcedPIshift);
end

if conjugate
    avg_output = conjugate(avg_output);
end

ofiles_mat = strrep(ofiles, 'tif', 'mat');
if ~five_bands
    Channels = norders;            % Depth = 3;
    avg_output_mat = avg_output(:, :, [1, 2:2:end]);
    save(ofiles_mat, 'avg_output_mat');
else
    Channels = norder * 2 - 1;     % Depth = 5;    
    save(ofiles_mat, 'avg_output');
end

header_output = ImageJ_formatted_TIFF.write_IFD(ostream_no, nz, fix(nx / 2) + 1, class(avg_output), Channels, 1, 1, dz, 0, 1, dr);
ImageJ_formatted_TIFF.write_IFD(ostream_exp_no, nz, fix(nx / 2) + 1, class(avg_output), Channels, 1, 1, dz, 0, 1, dr);

for i = 0 : (norders - 1)
    if i == 0
        Stack_band0 = permute(abs(avg_output(:, :, i + 1)), [2 1]);
        Stack_band0 = reshape(Stack_band0, [1, nz * (fix(nx / 2) + 1)]);
        Stack_band0_exp = Stack_band0;
        Stack_band0_exp(Stack_band0_exp >= 0.1) = 0.1;
        Stack_band0_exp = Stack_band0_exp ./ 0.1;
        Stack_band0_exp = Stack_band0_exp .^ 0.3;

        if system_endian ~= header_output.endian
            Stack_band0 = swapbytes(Stack_band0);
            Stack_band0_exp = swapbytes(Stack_band0_exp);
        end
        Stack_band0 = typecast(Stack_band0, 'uint8');
        Stack_band0_exp = typecast(Stack_band0_exp, 'uint8');
        fwrite(ostream_no, Stack_band0, 'uint8');
        fwrite(ostream_exp_no, Stack_band0_exp, 'uint8');
    else
        Stack_band0 = permute(abs(avg_output(:, :, 2 * i)), [2 1]);
        Stack_band0 = reshape(Stack_band0, [1, nz * (fix(nx / 2) + 1)]);
        Stack_band0_exp = Stack_band0;
        Stack_band0_exp(Stack_band0_exp >= 0.1) = 0.1;
        Stack_band0_exp = Stack_band0_exp ./ 0.1;
        Stack_band0_exp = Stack_band0_exp .^ 0.3;
        
        if system_endian ~= header_output.endian
            Stack_band0 = swapbytes(Stack_band0);
            Stack_band0_exp = swapbytes(Stack_band0_exp);
        end
        Stack_band0 = typecast(Stack_band0, 'uint8');
        Stack_band0_exp = typecast(Stack_band0_exp, 'uint8');
        fwrite(ostream_no, Stack_band0, 'uint8');
        fwrite(ostream_exp_no, Stack_band0_exp, 'uint8');
        
        if five_bands
            Stack_band0 = permute(abs(avg_output(:, :, 2 * i + 1)), [2 1]);
            Stack_band0 = reshape(Stack_band0, [1, nz * (fix(nx / 2) + 1)]);
            if system_endian ~= header_output.endian
                Stack_band0 = swapbytes(Stack_band0);
            end
            Stack_band0 = typecast(Stack_band0, 'uint8');
            fwrite(ostream_no, Stack_band0, 'uint8');
        end
    end
end
fclose(ostream_no);
fclose(ostream_exp_no);

if I2M_inc
    I2Mavg_output = fixorigin(I2Mavg_output, nx, nz, interpkr(1), interpkr(2));
    [I2Mavg_output, scalefactor] = rescale(I2Mavg_output, 0, nx, nz, scalefactor, dorescale);
    
    fileID = fopen(I2Mfiles, 'w+');
    header_I2Mfiles = ImageJ_formatted_TIFF.write_IFD(fileID, nz, fix(nx / 2 + 1), class(I2Mavg_output), 1, 1, 1, dz, 0, 1, dr);
    Stack_band0 = permute(abs(floatimage(:, :)), [2 1]);
    Stack_band0 = reshape(Stack_band0, [1, nz * (fix(nx / 2) + 1)]);
    
    if system_endian ~= header_I2Mfiles.endian
        Stack_band0 = swapbytes(Stack_band0);
    end
    Stack_band0 = typecast(Stack_band0, 'uint8');
    fwrite(fileID, Stack_band0, 'uint8');
    fclose(fileID);
end














%% Generates the matrix that is to be used to separate the raw indata into the different bands of sample information.
function sepMatrix = makematrix(nphases)
% INPUT VARIABLES
%   nphases: number of phases per direction, default = 5
% OUTPUT VARIABLE
%   	 phase = 0;         2pi / 5					2 * 2pi / 5				3 * 2pi / 5				4 * 2pi / 5
% sepMatrix = [1/5,			1/5,						1/5,						1/5,					1/5;				order = 0
%			  1/5,	cos(1 * 1 * 2pi / 5)/5,   cos(2 * 1 * 2pi / 5)/5,   cos(3 * 1 * 2pi / 5)/5,	cos(4 * 1 * 2pi / 5)/5;			order = 1, bandre
%			  0,	sin(1 * 1 * 2pi / 5)/5,   sin(2 * 1 * 2pi / 5)/5,   sin(3 * 1 * 2pi / 5)/5,	sin(4 * 1 * 2pi / 5)/5;			order = 1, bandim
%			  1/5,	cos(1 * 2 * 2pi / 5)/5,   cos(2 * 2 * 2pi / 5)/5,   cos(3 * 2 * 2pi / 5)/5,	cos(4 * 2 * 2pi / 5)/5;			order = 2, bandre
%			  0,	sin(1 * 2 * 2pi / 5)/5,   sin(2 * 2 * 2pi / 5)/5,   sin(3 * 2 * 2pi / 5)/5,	sin(4 * 2 * 2pi / 5)/5;]		order = 2, bandim

sepMatrix = zeros(nphases, nphases, 'single');

norders = fix((nphases + 1) / 2);											% norders = (5 + 1) / 2 = 3;
phi = 2 * pi / nphases;														% phi = 2pi / 5;
for j = 0 : (nphases - 1)
    sepMatrix(1, j + 1) = 1.0 / nphases;									% sepMatrix(1, 1:5) = 1.0 / 5;	order = 0;
    for order = 1 : (norders - 1)
        sepMatrix(2 * order, j + 1) = cos(j * order * phi) / nphases;		% sepMatrix(2, :), (4, :);
        sepMatrix(2 * order + 1, j + 1) = sin(j * order * phi) / nphases;	% sepMatrix(3, :), (5, :);
        % sepMatrix(2, :)(3, :) order = 1;
        % sepMatrix(4, :)(5, :) order = 2;
    end
end

end

%% Get correction files
function [background, slope] = getbg_and_slope(corrfiles, nx, ny)
% INPUT VARIABLES
%   corrfiles: correction file name
%   nx: 256
%   ny: 256
% OUTPUT VARIABLE
%   background: background image
%   slope: self-correction value for each pixel

[cstream_no, errmsg] = fopen(corrfiles, 'r');
if not(isempty(errmsg))
    warning(errmsg);
    error('File: %s does not exist.', ifiles);
end

header = ImageJ_formatted_TIFF.parse_tif(cstream_no, 0);
if (header.ImageWidth ~= nx) || (header.ImageLength ~= ny)
    error("Calibration file %s has different dimension than data file", corrfiles);
end

% calculate data_type from header.BitsPerSample
data_type = [header.SampleFormat, header.BitsPerSample];
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
    error("ReadTifStack does not support SampleFormat:%d with BitsPerSample: %d", header.SampleFormat, header.BitsPerSample);
end

% Initialize reading buffer and parameters
PixelNum = header.ImageWidth * header.ImageLength;
ByteCounts = fix(PixelNum * header.BitsPerSample / 8);
fseek(cstream_no, header.StripOffsets, 'bof');
background = typecast(transpose(fread(cstream_no, ByteCounts, 'uint8=>uint8')), dtype);
[~, ~, system_endian] = computer;
if system_endian ~= header.endian
    background = swapbytes(background);
end

background = reshape(background, [header.ImageWidth, header.ImageLength]);
background = single(permute(background, [2 1]));

% Calculate nz from header.ImageDescription and data_type from header.BitsPerSample
ImageDescription = convertStringsToChars(header.ImageDescription);
k1 = strfind(ImageDescription, 'images=');
ImageDescription_crop = ImageDescription(k1:end);
k2 = strfind(ImageDescription_crop, newline);
if ~isempty(k1) && ~isempty(k2)
    nz = uint32(str2double(ImageDescription(k1(1) + 7: k1(1) + k2(1) - 2)));
else
    warning("Did not find 'images=' in ImageDescription. Try to calculate it from filesize.");
    FileAttributes = dir(ifiles);
    nz = round(FileAttributes.bytes / (header.ImageWidth * header.ImageLength * (header.BitsPerSample / 8)));
end

% Correction file with only 1 section is assumed to contain background image
if nz > 1
    slope = typecast(transpose(fread(cstream_no, ByteCounts, 'uint8=>uint8')), dtype);
    if system_endian ~= header.endian
        slope = swapbytes(slope);
    end
    slope = reshape(slope, [header.ImageWidth, header.ImageLength]);
    slope = single(permute(slope, [2 1]));
else
    slope = ones([ny, nx], 'single');
end

fclose(cstream_no);
end

%% Estimate the background using the average intensities of the marginal pixels with border_size
function background = estimate_background(image, nx, ny, border_size)
% INPUT VARIABLES
%   image: single-layer input image
%   nx = header.nx = 256;
%   ny = header.ny = 256;
%   border_size = 20 pixels, default;
% OUTPUT VARIABLE
%   background


% C format
% total = 0;
% sum = 0;
%
% for i = 1: ny
%     if (i <= border_size || i > ny - border_size + 1)
%         for j = 1: nx
%             if (j <= border_size || j > nx - border_size + 1)
%                 sum = sum + image(i, j);
%                 total = total + 1;
%             end
%         end
%     end
% end
% background = single(sum / total);

% Matlab format
image_mask = ones(size(image), 'single');
image_mask(border_size + 1: ny - border_size, border_size + 1: nx - border_size) = 0;
background = sum(image .* image_mask, 'all') / sum(image_mask, 'all');

end

%% Softens the edges of a singe xy section to reduce edge artifacts and improve the fits
function image_apodized = apodize(napodize, nx, ny, image)
% INPUT VARIABLES
%   napodize = 10: number of pixels to soften the edges
%   nx = header.nx = 256;
%   ny = header.ny = 256;
%   image: input image;
% OUTPUT VARIABLE
%   image_apodized

image_apodized = image;

% Row-wise soften
diff = single(image_apodized(ny, :) - image_apodized(1, :)) / 2;    % Difference between last (256) row and the first (1) row;
for i = 1: napodize
    fact = single(1 - sin(((i - 0.5) / napodize) * pi * 0.5));
    image_apodized(i, :) = image_apodized(i, :) + diff .* fact;
    image_apodized(ny + 1 - i, :) = image_apodized(ny + 1 - i, :) - diff .* fact;
end

% Column-wise soften
diff = single(image_apodized(:, nx) - image_apodized(:, 1)) / 2;   % Difference between last (256) column and the first (1) column;
for j = 1: napodize
    fact = single(1 - sin(((j - 0.5) / napodize) * pi * 0.5));
    image_apodized(:, j) = image_apodized(:, j) + diff .* fact;
    image_apodized(:, nx + 1 - j) = image_apodized(:, nx + 1 - j) - diff .* fact;
end

image_apodized(image_apodized < 0) = 0;

end

%% Softens the edges of a singe xy section to reduce edge artifacts and improve the fits
function image_apodized = cosapodize(nx, ny, image)
% INPUT VARIABLES
%   nx = header.nx = 256;
%   ny = header.ny = 256;
%   image: input image;
% OUTPUT VARIABLE
%   image_apodized

% C format
% image_apodized = image;
% for j = 1:nx
%     xfact = single(sin(pi * (j - 0.5) / nx));
%     for i = 1:ny
%         yfact = single(sin(pi * (i - 0.5) / ny));
%         image_apodized(i, j) = image_apodized(i, j) * xfact * yfact;
%     end
% end

% Matlab format
x = 1:nx;
y = 1:ny;
[xfact, yfact] = meshgrid(x, y);
xfact = single(sin(pi * (xfact - 0.5) / nx));
yfact = single(sin(pi * (yfact - 0.5) / ny));
image_apodized = image .* xfact .* yfact;

image_apodized(image_apodized < 0) = 0;

end

%% Locate peak pixel to subpixel accuracy by fitting parabolas
function [xc, yc, zc, background, background_i2m, xc_i2m, yc_i2m, zc_i2m] = determine_center_and_background(stack5phases, I2M_image, nx, ny, nz, nphases, twolens, I2M)
% INPUT VARIABLES
%   stack5phases: input image
%   nx = header.nx = 256;
%   ny = header.ny = 256;
%   nphases = 5;
%   twolens = 0;    Whether using I5S two objective acquisition
% OUTPUT VARIABLE
%   xc: 
%   yc: 
%   zc: 
%   background:

fprintf("Locate PSF center using parabolic fitting.\n");
% stack3d = zeros(ny, nx ,nz, 'single');

% stack3d now is the average of 5 phases images;
stack3d = single(sum(stack5phases, 4) ./ nphases);

% Whether to apply 3-D Gaussian filtering
if (false)
    stack3d = imgaussfilt3(stack3d, 1);
    fileID = fopen('stack_filtered.tif', 'w+');
    [Height, Width, Depth] = size(stack3d);
    minval = min(stack3d, [], 'all');
    maxval = max(stack3d, [], 'all');
    header_bands = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(stack3d), 1, Depth, 1, 0.06, minval, maxval, 0.0709);
    Stack_bands = permute(stack3d, [2 1 3]);
    Stack_bands = reshape(Stack_bands, [1, nx * ny * nz]);
    
    [~, ~, system_endian] = computer;
    if system_endian ~= header_bands.endian
        Stack_bands = swapbytes(Stack_bands);
    end
    Stack_bands = typecast(Stack_bands, 'uint8');
    fwrite(fileID, Stack_bands, 'uint8');
    fclose(fileID);
end

% Search for the peak pixel
% Be aware that stack3d is of dimension [ny, nx, nz]
[maxval, idx] = max(stack3d(:));
[maxi, maxj, maxk] = ind2sub(size(stack3d), idx);
disp(['Maximum pixel intensity: ', num2str(maxval)]);
disp(['y index: ', num2str(maxi)]);
disp(['x index: ', num2str(maxj)]);
disp(['z index: ', num2str(maxk)]);


iminus = maxi - 1; iplus = maxi + 1;
if (iminus <= 0) 
    iminus = iminus + ny; end
if (iplus > ny)
    iplus = iplus - ny; end

jminus = maxj - 1; jplus = maxj + 1;
if (jminus <= 0)
    jminus = jminus + nx; end
if (jplus > nx)
    jplus = jplus - nx; end

kminus = maxk - 1; kplus = maxk + 1;
if (kminus <= 0)
    kminus = kminus + nz; end
if (kplus > nz)
    kplus = kplus - nz; end

% find z center of mass
if ~twolens
        valminus = stack3d(maxi, maxj, kminus);
        valplus = stack3d(maxi, maxj, kplus);
        % find z center of mass using 3 points parabolic fitting
        zc = maxk + fitparabola(valminus, maxval, valplus);
else
    if (true)
        % find z center of mass using 3 points parabolic fitting
        zcs = zeros([5, 1], 'single'); peakval = zeros([5, 1], 'single');
        n_samples = 3;
        Xs = zeros([n_samples, 1], 'single'); Ys = zeros([n_samples, 1], 'single');
        
        % 3 points: kminus2 = maxk - 1; kplus2 = maxk + 1;
        kminus2 = maxk - fix(n_samples / 2); kplus2 = maxk + fix(n_samples / 2);
        
        % search for the current trough
        for i = 1:n_samples
            % 3 points: Xs(1:3) = maxk -1 +0 +1;
            Xs(i) = kminus2 + i - 1;
            Ys(i) = stack3d(maxi, maxj, Xs(i));
        end
        [zcs(1), peakval(1)] = fitparabola_Nsample_points(Xs, Ys, n_samples);
        
        minval_right = Ys(n_samples);       % minval_right = Ys(3) = stack3d(maxi, maxj, Xs(3)) = stack3d(maxi, maxj, maxk + 1)
        minval_left = Ys(1);                % minval_left = stack3d(maxi, maxj, Xs(1)) = stack3d(maxi, maxj, maxk - 1)
        
        % search forward for the next trough
        kmin = kplus2;                      % kmin = maxk + 1;
        if stack3d(maxi, maxj, kplus2 + 1) < minval_right
            minval_right = stack3d(maxi, maxj, kplus2 + 1);
            kmin = kplus2 + 1;              % kmin = maxk + 2;
        end

        for i = 1:n_samples
            Xs(i) = kmin + i - 1;           % Xs(1:3) = maxk +2 +3 +4;
            Ys(i) = stack3d(maxi, maxj, Xs(i));
        end
        [zcs(2), peakval(2)] = fitparabola_Nsample_points(Xs, Ys, n_samples);
        
        % search backward for the previous trough
        kmin = kminus2;                     % kmin = maxk - 1;
        if stack3d(maxi, maxj, kminus2 - 1) < minval_left
            minval_left = stack3d(maxi, maxj, kminus2 - 1);
            kmin = kminus2 - 1;             % kmin = maxk - 2;
        end
        for i = 1:n_samples
            Xs(i) = kmin - i + 1;           % Xs(1:3) = maxk -2 -3 -4;
            Ys(i) = stack3d(maxi, maxj, Xs(i));
        end
        [zcs(3), peakval(3)] = fitparabola_Nsample_points(Xs, Ys, n_samples);
        
        
    else
        % find z center of mass using 5 points parabolic fitting stead of 3 points
        zcs = zeros([5, 1], 'single'); peakval = zeros([5, 1], 'single');
        n_samples = 5;
        Xs = zeros([n_samples, 1], 'single'); Ys = zeros([n_samples, 1], 'single');
        
        % 5 points: kminus2 = maxk - 2; kplus2 = maxk + 2;
        kminus2 = maxk - fix(n_samples / 2); kplus2 = maxk + fix(n_samples / 2);
        
        % search for the current trough
        for i = 1:n_samples
            % 5 points: Xs(1:5) = maxk -2 -1 +0 +1 +2;
            Xs(i) = kminus2 + i - 1;
            Ys(i) = stack3d(maxi, maxj, Xs(i));
        end
        [zcs(1), peakval(1)] = fitparabola_Nsample_points(Xs, Ys, n_samples);
        
        minval_right = Ys(n_samples);       % minval_right = Ys(5) = stack3d(maxi, maxj, Xs(5)) = stack3d(maxi, maxj, maxk + 2)
        minval_left = Ys(1);                % minval_left = stack3d(maxi, maxj, Xs(1)) = stack3d(maxi, maxj, maxk - 2)
        
        % search forward for the next trough
        kmin = kplus2;                              % kmin = maxk + 2;
        for k = kplus2 + 1: kplus2 + n_samples      % k = maxk + 3: maxk + 7
            if stack3d(maxi, maxj, k) < minval_right
                minval_right = stack3d(maxi, maxj, k);
                kmin = k;                   % kmin = (maxk + 3: maxk + 7), suppose kmin = maxk + 5;
            end
        end
        for i = 1:n_samples
            Xs(i) = kmin + i;               % Xs(1:5) = maxk +6 +7 +8 +9 +10;
            Ys(i) = stack3d(maxi, maxj, Xs(i));
        end
        [zcs(2), peakval(2)] = fitparabola_Nsample_points(Xs, Ys, n_samples);
        
        % search backward for the previous trough
        kmin = kminus2;                             % kmin = maxk - 2;
        for k = kminus2 - 1:-1:kminus2 - n_samples  % k = maxk - 7: maxk - 3
            if stack3d(maxi, maxj, k) < minval_left
                minval_left = stack3d(maxi, maxj, k);
                kmin = k;                   % kmin = (maxk - 7: maxk - 3), suppose kmin = maxk - 5;
            end
        end
        for i = 1:n_samples
            Xs(i) = kmin - i;               % Xs(1:5) = maxk -6 -7 -8 -9 -10;
            Ys(i) = stack3d(maxi, maxj, Xs(i));
        end
        [zcs(3), peakval(3)] = fitparabola_Nsample_points(Xs, Ys, n_samples);
    end
    
    fprintf("zc0=%.3f, zc_plus=%.3f, zc_minus=%.3f\n", zcs(1), zcs(2), zcs(3));
    
    % find z center of mass using 3 points parabolic fitting
    [zc, ~] = fitparabola_Nsample_points(zcs(1:3), peakval(1:3), 3);
end

% find y center of mass
valminus = stack3d(iminus, maxj, maxk);
valplus = stack3d(iplus, maxj, maxk);
yc = maxi + fitparabola(valminus, maxval, valplus);

% find x center of mass
valminus = stack3d(maxi, jminus, maxk);
valplus = stack3d(maxi, jplus, maxk);
xc = maxj + fitparabola(valminus, maxval, valplus);

% Caluclate background
infocus_sec = floor(zc);
background = sum(stack5phases(1:floor(yc) - 20, :, infocus_sec, :), 'all') / (nphases * (floor(yc) - 20) * nx);

if I2M
    background_i2m = sum(I2M_image(1:floor(yc) - 20, :, infocus_sec, :), 'all') / (nphases * (floor(yc) - 20) * nx);
    [maxval, idx] = max(I2M_image(:));
    [maxi, maxj, maxk] = ind2sub(size(I2M_image), idx);
    iminus = maxi - 1; iplus = maxi + 1;
    
    if (iminus <= 0)
        iminus = iminus + ny; end
    if (iplus > ny)
        iplus = iplus - ny; end
    
    jminus = maxj - 1; jplus = maxj + 1;
    if (jminus <= 0)
        jminus = jminus + nx; end
    if (jplus > nx)
        jplus = jplus - nx; end
    
    kminus = maxk - 1; kplus = maxk + 1;
    if (kminus <= 0)
        kminus = kminus + nz; end
    if (kplus > nz)
        kplus = kplus - nz; end
    
    % find z center of mass
    valminus = I2M_image(maxi, maxj, kminus);
    valplus = I2M_image(maxi, maxj, kplus);
    zc_i2m = maxk + fitparabola(valminus, maxval, valplus);
    
    % find y center of mass
    valminus = I2M_image(iminus, maxj, maxk);
    valplus = I2M_image(iplus, maxj, maxk);
    yc_i2m = maxi + fitparabola(valminus, maxval, valplus);
    
    % find x center of mass
    valminus = I2M_image(maxi, jminus, maxk);
    valplus = I2M_image(maxi, jplus, maxk);
    xc_i2m = maxj + fitparabola(valminus, maxval, valplus);
else
    background_i2m = [];
    xc_i2m = []; yc_i2m = []; zc_i2m = [];
end
end

%% Locate peak pixel to subpixel accuracy by using window binning maximum
function [xc, yc, zc, background, background_i2m, xc_i2m, yc_i2m, zc_i2m] = determine_center_and_background_2(stack5phases, I2M_image, nx, ny, nz, nphases, I2M)
% INPUT VARIABLES
%   stack5phases: input image
%   nx = header.nx = 256;
%   ny = header.ny = 256;
%   nphases = 5;
% OUTPUT VARIABLE
%   xc: 
%   yc: 
%   zc: 
%   background:

fprintf("Locate PSF center using window binning maximum.\n");
% stack3d = zeros(ny, nx ,nz, 'single');

% stack3d now is the average of 5 phases images;
stack3d = single(sum(stack5phases, 4) ./ nphases);


win1 = 2;
win2 = 5;
win21 = win1 + win2;

% locate peak intensity in (2 * win1 + 1)^3 window
sum_matrix = zeros([ny - 2 * win21, nx - 2 * win21, nz - 2 * win21], 'single');

for dz = -win1 : win1
    for dy = -win1 : win1
        for dx = -win1 : win1
            sum_matrix = sum_matrix + stack3d(1 + win21 + dy : dy - win21 + ny, 1 + win21 + dx : dx - win21 + nx, 1 + win21 + dz : dz - win21 + nz);
        end
    end
end

% Search for the peak pixel
% Be aware that sum_matrix is of dimension [ny - 2 * win21, nx - 2 * win21, nz - 2 * win21]
[maxval, idx] = max(sum_matrix(:));
[maxi, maxj, maxk] = ind2sub(size(sum_matrix), idx);
maxi = maxi + win21;
maxj = maxj + win21;
maxk = maxk + win21;
disp(['Maximum pixel intensity after 5 by 5 by 5 binning: ', num2str(maxval)]);
disp(['y index: ', num2str(maxi)]);
disp(['x index: ', num2str(maxj)]);
disp(['z index: ', num2str(maxk)]);

% do a center-of-mass fit in (2 * win2 + 1)^3 window
pos = zeros([4, 1], 'single');

for y = maxi - win2  : maxi + win2
    for x = maxj - win2 : maxj + win2
        for z = maxk - win2 : maxk + win2
            value = stack3d(y, x, z);
            pos(1) = pos(1) + value;
            pos(2) = pos(2) + value * x;
            pos(3) = pos(3) + value * y;
            pos(4) = pos(4) + value * z;
        end
    end
end

xc = pos(2) / pos(1);
yc = pos(3) / pos(1);
zc = pos(4) / pos(1);

% Caluclate background
infocus_sec = floor(zc);
background = sum(stack5phases(1:floor(yc) - 20, :, infocus_sec, :), 'all') / (nphases * (floor(yc) - 20) * nx);

if I2M
    background_i2m = sum(I2M_image(1:floor(yc) - 20, :, infocus_sec, :), 'all') / (nphases * (floor(yc) - 20) * nx);
    % locate peak intensity in win1^3 window
    sum_matrix = zeros([ny - 2 * win21, nx - 2 * win21, nz - 2 * win21], 'single');
    
    for dz = -win1 : win1
        for dy = -win1 : win1
            for dx = -win1 : win1
                sum_matrix = sum_matrix + I2M_image(1 + win21 + dy : dy - win21 + ny, 1 + win21 + dx : dx - win21 + nx, 1 + win21 + dz : dz - win21 + nz);
            end
        end
    end
    
    % Search for the peak pixel
    % Be aware that sum_matrix is of dimension [ny - 2 * win21, nx - 2 * win21, nz - 2 * win21]
    [maxval, idx] = max(sum_matrix(:));
    [maxi, maxj, maxk] = ind2sub(size(sum_matrix), idx);
    maxi = maxi + win21;
    maxj = maxj + win21;
    maxk = maxk + win21;
    disp(['Maximum pixel intensity: ', num2str(maxval)]);
    disp(['y index: ', num2str(maxi)]);
    disp(['x index: ', num2str(maxj)]);
    disp(['z index: ', num2str(maxk)]);
    
    % do a center-of-mass fit in a window +-4
    pos = zeros([4, 1], 'single');
    
    for y = maxi - win2  : maxi + win2
        for x = maxj - win2 : maxj + win2
            for z = maxk - win2 : maxk + win2
                value = I2M_image(y, x, z);
                pos(1) = pos(1) + value;
                pos(2) = pos(2) + value * x;
                pos(3) = pos(3) + value * y;
                pos(4) = pos(4) + value * z;
            end
        end
    end
    
    xc_i2m = pos(2) / pos(1);
    yc_i2m = pos(3) / pos(1);
    zc_i2m = pos(4) / pos(1);
else
    background_i2m = [];
    xc_i2m = []; yc_i2m = []; zc_i2m = [];
end


end

%% Fits a parabola to the three points (-1,a1), (0,a2), and (1,a3).
function peak = fitparabola(a1, a2, a3)
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

%% Fits a parabola to the N points (Xs(1), Ys(1)), (Xs(2), Ys(2), ..., (Xs(N), Ys(N))
function [peak, value] = fitparabola_Nsample_points(Xs, Ys, n)
% INPUT VARIABLES
%   Xs:
%   Ys:
%   n:
% OUTPUT VARIABLE
%   peak: x-value of the max (or min) of the parabola
%   value: y-value of the max (or min) of the parabola

p = polyfit(Xs, Ys, 2);
% p(x) = p(1)*x^2 + p(2)*x + p(3)
peak = single(-p(2) / (2 * p(1)));
value = single(p(3) - p(2)^2 / (4 * p(1)));
end

%% Separate floatimage into real and imaginary bands
function output = separate(nx, ny, z, nphases, floatimage, sepMatrix)
% Applies image arithmetic to the image sequence floatimage[][],
% which was acquired with different phases of the illumination
% pattern, to separate the different bands of sample information.
% The coefficients are pre-stored in the matrix sepMatrix.
% The bands are returned in the same array floatimage where the
% input data was.
% INPUT VARIABLES
%   nx = header.nx = 256;
%   ny = header.ny = 256;
%   z: current z index;
% OUTPUT VARIABLE
%   output

MAXPHASES = 10;
if (nphases > MAXPHASES)
    error("In separate(), nphases is larger than MAXPHASES\n");
end

output = zeros(ny, nx, nphases, 'single');

% Matlab format
for i = 1 : nphases
    for j = 1 : nphases
        output(:, :, i) = output(:, :, i) + floatimage(:, :, z, j) .* sepMatrix(i, j);
    end
end
% output(:, :, order = 0)         = Sum(sepMatrix(order = 0, phase= 1,2,3,4,5) * floatimage(:, :, z, phase = 1,2,3,4,5);
% output(:, :, order = 1, bandre) = Sum(sepMatrix(order = 1, bandre, phase= 1,2,3,4,5) * floatimage(:, :, z, phase = 1,2,3,4,5);
% output(:, :, order = 1, bandim) = Sum(sepMatrix(order = 1, bandim, phase= 1,2,3,4,5) * floatimage(:, :, z, phase = 1,2,3,4,5);
% output(:, :, order = 2, bandre) = Sum(sepMatrix(order = 2, bandre, phase= 1,2,3,4,5) * floatimage(:, :, z, phase = 1,2,3,4,5);
% output(:, :, order = 2, bandim) = Sum(sepMatrix(order = 2, bandim, phase= 1,2,3,4,5) * floatimage(:, :, z, phase = 1,2,3,4,5);

% output(:, :, order = 0) = 1/5 * Sum(floatimage(:, :, :, phase = 1, 2, 3, 4, 5));
% output(:, :, order = 1, bandre) = [1/5, cos(1 * 1 * 2pi / 5)/5, cos(2 * 1 * 2pi / 5)/5, cos(3 * 1 * 2pi / 5)/5, cos(4 * 1 * 2pi / 5)/5] *
%                                   [floatimage(1); floatimage(2); floatimage(3); floatimage(4); floatimage(5)]
% output(:, :, order = 1, bandim) = ...
% output(:, :, order = 2, bandre) = ...
% output(:, :, order = 2, bandim) = ...
end

%% Modify the phase of bands, so that it corresponds to FFT of a bead at origin
function bands_return = shift_center(bands, nx, ny, nz, xc, yc, zc)
% To get rid of checkerboard effect in the OTF bands
% (xc, yc, zc) is the estimated center of the point source, which in most cases is the bead
% Converted from Fortran code. kz is treated differently than kx and ky. don't know why
% INPUT VARIABLES
%   bands:
%   nx = header.nx = 256;
%   ny = header.ny = 256;
%   nz = header.nz / nphases = 64;
%   xc: supposed as 127.2
%   yc: supposed as 128.5
%   zc: supposed as 32.3
% OUTPUT VARIABLE
%   bands_return:


kycent = fix(ny / 2);                   % kycent = 128;
kxcent = fix(nx / 2);                   % kxcent = 128;
kzcent = fix(nz / 2);                   % kzcent = 32;

% the origin of spatial domain is at (1, 1, 1), and we want to move PSF center
% from (yc, xc, zc) to (1, 1, 1)
dphiz = 2 * pi * (zc - 1) / nz;               % dphiz = 2 * pi * 31.3 / 64;
dphiy = 2 * pi * (yc - 1) / ny;               % dphiy = 2 * pi * 127.5 / 256;
dphix = 2 * pi * (xc - 1) / nx;               % dphix = 2 * pi * 126.2 / 256;

% C format
% bands_return = complex(zeros([ny, nx, nz], 'single'));
% for kin = 0:(nz - 1)
%     kz = kin;
%     if kz > kzcent
%         kz = kz - nz;
%     end
%     phi1 = dphiz * kz;                  % first part of phi (z offset)
%         for iin = 0:(ny - 1)
%             ky = iin;
%             if iin > kycent
%                 ky = ky - ny;
%             end
%             phi2 = dphiy * ky;          % second part of phi (y offset)   
%             for jin = 0:(nx - 1)
%                 kx = jin;
%                 if kx > kxcent
%                     kx = kx - nx;
%                 end
%                 phi3 = dphix * kx;      % third part of phi (x offset)
%                 phi = phi1 + phi2 + phi3;                
%                 bands_return(iin + 1, jin + 1, kin + 1) = single(bands(iin, jin, kin) * exp(1i * phi));
%             end
%         end
% end

% Matlab format
y = linspace(0, ny - 1, ny);
y(y > kycent) = y(y > kycent) - ny;
x = linspace(0, nx - 1, nx);
x(x > kxcent) = x(x > kxcent) - nx;
z = linspace(0, nz - 1, nz);
z(z > kzcent) = z(z > kzcent) - nz;

[X, Y, Z] = meshgrid(x, y, z);

phi1 = dphiz .* Z;
phi2 = dphiy .* Y;
phi3 = dphix .* X;
phi = phi1 + phi2 + phi3;
bands_return = single(bands .* exp(1i * phi));

end

%% Dividing the measured OTF, before rotational averaging, by S(k - mp)
function bands_return = beadsize_compensate(bands, k0angle, linespacing, bead_diameter, norders, nx, ny, nz, dkr, dkz)
% To compensate finite bead size for all bands
% INPUT VARIABLES
%   bands:
%   k0angle = k0angleguess = 0.008 rad; It should be noted that this
%   k0angle is opposite to the Cartesian coordinate system in C, matlab (y direction)
%   "inverse" to that we 
%   linespacing = 0.198802 um;
%   bead_diameter = 0.109 um;
%   norders = (nphases + 1) / 2 = 3;
%   nx = header.nx = 256;
%   ny = header.ny = 256;
%   nz = header.nz / nphases = 64;
%   dkr = 1 / (ny * dr) = 1 / (256 * 0.08) = 1 / 20.48 um;
%   dkz = 1 / (nz * dz) = 1 / (64 * 0.125) = 1 / 8 um;
% OUTPUT VARIABLE
%   bands_return:

fprintf("In beadsize_compensate()\n");
bands_return = zeros(size(bands), 'single');
kycent = fix(ny / 2);                   % kycent = 128;
kxcent = fix(nx / 2);                   % kxcent = 128;
kzcent = fix(nz / 2);                   % kzcent = 32;

k0mag = single(1 / linespacing);                % k0mag = 1 / 0.198802 um;
radius = single(bead_diameter * 0.5);           % radius = 0.109 um * 0.5 = 0.0545 um;
% the limit value of the FT of a sphere at the origin of Fourier space
limit_at_origin = single(4 * pi * radius ^ 3 / 3);    % limit_at_origin = 4 / 3 * pi * (0.0545 um) ^ 3;

for order = 0 : (norders - 1)
    if order == 0
        bandptr = bands(:, :, :, order + 1);        % order = 0
    else
        bandptr = bands(:, :, :, 2 * order);        % order = 1 or 2, bandre
        bandptr1 = bands(:, :, :, 2 * order + 1);   % order = 1 or 2, bandim
    end    
    
    %{
    (k0x, k0y) are the frequency components along x, y axis
    When order = 0: (k0x, k0y) = (0, 0);
    When order = 1: (k0x, k0y) = (1 / (2 * 0.198802 um) * cos(0.008), 1 / (2 * 0.198802 um) * sin(0.008));
    When order = 2: (k0x, k0y) = (2 / (2 * 0.198802 um) * cos(0.008), 2 / 2 * 0.198802 um * sin(1.57193));
    %}
    k0x = single(order / (norders - 1) * k0mag * cos(k0angle));
    k0y = single(order / (norders - 1) * k0mag * sin(k0angle));
    
    fprintf("order=%d, k0x=%f, k0y=%f\n", order, k0x, k0y);
    
    % C format
%     for kin = 0 : (nz - 1)          % kin = 0:63
%         kout = kin;
%         if kout > kzcent
%             kout = kout - nz;       % kout = -31:32
%         end
%         kz = kout * dkz;            % kz = (-31:32) / 8 um        
%         for iin = 0 : (ny - 1)      % iin = 0:255
%             iout = iin;
%             if iout > kycent
%                 iout = iout - ny;   % iout = -127:128
%             end
%             ky = iout * dkr + k0y;  % ky = (-127:128) / 20.48 um + k0y
%             for jin = 0 : (nx - 1)
%                 jout = jin;
%                 if jout > kxcent
%                     jout = jout - nx;
%                 end
%                 kx = jout * dkr + k0x;  % kx = (-127:128) / 20.48 um + k0x
%                 
%                 if ~((order == 0) && (kin == 0) && (iin == 0) && (jin == 0))
%                     % rho: the distance in Fourier space from a point to the origin;
%                     rho = sqrt(kz * kz + ky * ky + kx * kx);
%                     ratio = sphereFFT(rho, radius) / limit_at_origin;
%                 else
%                     ratio = 1;
%                 end
%                 
%                 if order == 0
%                     bandptr(iin + 1, jin + 1, kin + 1) = bandptr(iin + 1, jin + 1, kin + 1) ./ ratio;
%                     bands_return(iin + 1, jin + 1, kin + 1, order + 1) = bandptr(iin + 1, jin + 1, kin + 1);
%                 else
%                     bandptr(iin + 1, jin + 1, kin + 1) = bandptr(iin + 1, jin + 1, kin + 1) ./ ratio;
%                     bandptr1(iin + 1, jin + 1, kin + 1) = bandptr1(iin + 1, jin + 1, kin + 1) ./ ratio;
%                     bands_return(:, :, :, 2 * order) = bandptr(iin + 1, jin + 1, kin + 1);
%                     bands_return(:, :, :, 2 * order + 1) = bandptr1(iin + 1, jin + 1, kin + 1);
%                 end                
%             end
%         end
%     end
    
    % Matlab format
    y = linspace(0, ny - 1, ny);
    y(y > kycent) = y(y > kycent) - ny;
    ky = single(y .* dkr + k0y);
    x = linspace(0, nx - 1, nx);
    x(x > kxcent) = x(x > kxcent) - nx;
    kx = single(x .* dkr + k0x);
    z = linspace(0, nz - 1, nz);
    z(z > kzcent) = z(z > kzcent) - nz;
    kz = single(z .* dkz);
    
    [kX, kY, kZ] = meshgrid(kx, ky, kz);
    rho = sqrt(kZ .* kZ + kY .* kY + kX .* kX);
    ratio = sphereFFT(rho, radius) / limit_at_origin;   % When rho = 0, ratio = 1.
    
    if order == 0
        bandptr = bandptr ./ ratio;
        bands_return(:, :, :, order + 1) = bandptr;
    else
        bandptr = bandptr ./ ratio;
        bandptr1 = bandptr1 ./ ratio;
        bands_return(:, :, :, 2 * order) = bandptr;
        bands_return(:, :, :, 2 * order + 1) = bandptr1;
    end    
end

end

%% Fourier transform of a solid sphere with a diameter equal to the nominal diameter of the bead
function a = sphereFFT(k, radius)
% According to Mats's derivation,
% FT of a sphere with radius R is f(k) = (R / Pi * k ^ 2) * g(2 * Pi * k * R),
% where g(x) = sin(x) / x - cos(x)
% The limit at the origin of F space is 4 * Pi * R ^ 3 / 3

% C format
% if (k > 0 || k < 0)
%     x = 2 * pi * radius * k;
%     a = radius / (pi * k * k) * (sin(x) / x - cos(x));
% else
%     a = 4 * pi * (radius ^ 3) / 3;
% end

% Matlab format
a = zeros(size(k), 'single');
x = k;
x(k ~= 0) = 2 * pi * radius .* k(k ~= 0);
a(k ~= 0) = radius ./ (pi .* k(k ~= 0) .^ 2) .* (sin(x(k ~= 0)) ./ x(k ~= 0) - cos(x(k ~= 0)));
a(k == 0) = 4 * pi * (radius ^ 3) / 3;

end

%% To decrease noise, each Om was rotationally averaged around the kz axis
function avg_output = radialft(band, nx, ny, nz)
% To decrease noise, each Om was rotationally averaged around the kz axis
% INPUT VARIABLES
%   bands:
%   nx = header.nx = 256;
%   ny = header.ny = 256;
%   nz = header.nz / nphases = 64;
% OUTPUT VARIABLE
%   avg_output:

fprintf("In radialft()\n");
kycent = fix(ny / 2);                   % kycent = 128;
kxcent = fix(nx / 2);                   % kxcent = 128;
kzcent = fix(nz / 2);                   % kzcent = 32;
avg_output = complex(zeros(kxcent + 1, nz, 'single')); % [129, 64]
count = zeros(kxcent + 1, nz, 'single');               % [129, 64]

% Accumulate values at radius rdist (xy plane)
% C format
% for kin = 0 : (nz - 1)
%     if kin > kzcent
%         kz = kz - nz;
%     end
%     for iin = 0 : (ny - 1)
%         ky = iin;
%         if iin > kycent
%             ky = ky - ny;
%         end
%         for jin = 0 : kxcent
%             kx = jin;
%             rdist = sqrt(kx ^ 2 + ky ^ 2);
%             if round(rdist) < kxcent + 1
%                 avg_output(round(rdist) + 1, kin + 1) = avg_output(round(rdist) + 1, kin) + band(iin + 1, jin + 1, kin + 1);
%                 count(round(rdist) + 1, kin + 1) = count(round(rdist) + 1, kin + 1) + 1;
%             end
%         end
%     end
% end

% Matlab format
for iin = 0 : (ny - 1)
    ky = iin;
    if ky > kycent
        ky = ky - ny;
    end    
    for jin = 0 : kxcent
        kx = jin;
        rdist = sqrt(kx ^ 2 + ky ^ 2);
        if round(rdist) < kxcent + 1
            avg_output(round(rdist) + 1, :) = avg_output(round(rdist) + 1, :)...
                + squeeze(band(iin + 1, jin + 1, :)).';
            count(round(rdist) + 1, :) = count(round(rdist) + 1, :) + 1;
        end
    end
end

% Divided by count
% C format
% for kin = 1 : nz
%     for jin = 1 : kxcent + 1
%         avg_output(jin, kin) = avg_output(jin, kin) / count(jin, kin);
%     end
% end

% Matlab format
avg_output =  avg_output ./ count;

% Then complete the rotational averaging and scaling in z axis
% C format
% for kx = 0 : kxcent
%     avg_output(kx + 1, 1) = real(avg_output(kx + 1, 1));    
%     for kz = 1 : kzcent
%         indout = kz + 1;
%         indout_conj = nz - kz + 1;
%         avg_output_real = (real(avg_output(kx + 1, indout)) + real(avg_output(kx + 1, indout_conj))) / 2;
%         avg_output_imag = (imag(avg_output(kx + 1, indout)) - imag(avg_output(kx + 1, indout_conj))) / 2;
%         avg_output(kx + 1, indout) = complex(avg_output_real, avg_output_imag);
%         avg_output(kx + 1, indout_conj) = conj(avg_output(kx + 1, indout));
%     end
% end

% Matlab format
avg_output(:, 1) = real(avg_output(:, 1));              % DC component is always real
for kz = 1 : kzcent
    indout = kz + 1;
    indout_conj = nz - kz + 1;
    avg_output_real = (real(avg_output(:, indout)) + real(avg_output(:, indout_conj))) / 2;
    avg_output_imag = (imag(avg_output(:, indout)) - imag(avg_output(:, indout_conj))) / 2;
    avg_output(:, indout) = complex(avg_output_real, avg_output_imag);
    avg_output(:, indout_conj) = conj(avg_output(:, indout));
%     avg_output(:, indout) = (avg_output(:, indout) + conj(avg_output(:, indout_conj))) / 2;
%     avg_output(:, indout_conj) = conj(avg_output(:, indout));
end

end

%% Modify the kx = 0 row with kx = 1 row; modify the kz = 0 column with average of kz = 2 and last (ny - 1) column (partially)
function avg_output = modify(otfkxkz, nx, nz, ifixz, ifixr, order, twolens)
% INPUT VARIABLES
%   otfkxkz = avg_output(:, :, phase);
%   nx = header.nx = 256;
%   nz = header.nz / nphases = 64;
%   ifixz = [1, 1, 1];
%   ifixr = 0;
%   order = fix(phase / 2) = 0, 1, 2;
%   twolens = 0;    Whether using I5S two objective acquisition
% OUTPUT VARIABLE
%   avg_output

avg_output = otfkxkz;
if ifixz(1)
    % C format
%     for kz = 0 : (nz - 1)
%         if (kz ~= 0) || (order > 0)
%             avg_output(1, kz + 1) = avg_output(2, kz + 1);
%         end
%     end
    
    % Matlab format
    if order == 0
        avg_output(1, 2 : end) = avg_output(2, 2 : end);
    elseif order > 0
        avg_output(1, :) = avg_output(2, :);
    end
    % replace the kx = 0 row with kx = 1 row
end

if ifixr > 0
    % C format
%     for kx = ifixr : fix(nx / 2)
%         avg_output(kx + 1, 1) = (avg_output(kx + 1, 2) + avg_output(kx + 1, nz)) / 2;
%     end
    
    % Matlab format
    kx = ifixr + 1 : fix(nx / 2) + 1;
    avg_output(kx, 1) = (avg_output(kx, 2) + avg_output(kx, nz)) / 2;
    % replace 1st column within range (ifixr + 1:end) with average of 2nd
    % column and last column    
end
end

%% Set to zero outside of OTF's known support
function avg_output = cleanup(otfkxkz, order, nx, nz, dkr, dkz, linespacing,...
    lamdanm, icleanup, ileave0_1, ileave0_2, ileave1_1, ileave1_2, ileave1_3, ileave1_4, ileave2, twolens, NA, NIMM)
% INPUT VARIABLES
%   otfkxkz = avg_output(:, :, phase);
%   order = fix(phase / 2) = 0, 1, 2;
%   nx = header.nx = 256;
%   nz = header.nz / nphases = 64;
%   dkr = 1 / (ny * dr) = 1 / (256 * 0.08) = 1 / 20.48 um;
%   dkz = 1 / (nz * dz) = 1 / (64 * 0.125) = 1 / 8 um;
%   linespacing = 0.198802 um;
%   icleanup = fix(nx / 2) + 1 = 129;
%   lamdanm = 525 nm
%   ileave1_1 = ileavekz(1) = 16; ileave1_1 is the minimum shifted-pixel number for order1 OTF;
%   ileave1_2 = ileavekz(2) = 20; ileave1_2 is the maximum shifted-pixel number for order1 OTF;
%   ileave2 = ileavekz(3) = 3;    ileave2 is the maximum axially-broadened pixel number for order2 OTF;
%   twolens = 0;    Whether using I5S two objective acquisition
%   NA = 1.35
%   NIMM = 1.406;
% OUTPUT VARIABLE
%   avg_output

SPOTRATIO = 0.1;
avg_output = otfkxkz;

NA_local = NA;                          % NA_local = NA = 1.35;
lamda = lamdanm * 0.001;                % lamda = 525 nm * 0.001 = 0.525 um;
sinalpha = NA_local / NIMM;             % sinalpha = 1.35 / 1.406;
cosalpha = cos(asin(sinalpha));         % cosalpha = cos(asin(1.35 / 1.406));
krmax = 2 * NA_local / lamda;           % krmax = 2 * 1.35 / (0.525 um) = 5.142857 um-1; The inverse form of Abbe limit;
                                        % The maximum lateral frequency of NA and lamda;
                                        % krmax / dkr = 105.3 pixels
k0mag = single(1.0 / linespacing);      % k0mag = 1.0 / (0.189 um);

fprintf("nz=%d\n", nz);

if ~twolens
    for ix = 0 : (icleanup - 1)         % ix = 0:128
        kr = ix * dkr;                  % kr = 0:6.25
        % ix = 0:105
        if kr <= krmax
            beta = asin((NA_local - kr * lamda) / NIMM);
            % When kr = 0:          beta = asin(NA_local / NIMM) = asin(1.35 / 1.406) = 1.2876 (73.7748 degrees);
            % When kr = krmax / 2:  beta = 0
            % When kr = krmax:      beta = -1.2876 (-73.7748 degrees);
            % It can be regarded as: beta = asin((krmax / 2 - kr) / (NIMM / lamda));
            kzedge = (NIMM / lamda) * (cos(beta) - cosalpha);
%             fprintf("kzedge=%f, dkz=%f\n", kzedge, dkz);
            if order == 0
                kzstart = round((kzedge / dkz) + 1);
                kzend = nz - kzstart;
                % When kr = 0 or krmax: kzedge = (NIMM / lamda) * (cosalpha - cosalpha) = 0
                %   kzstart = round(0 / dkz + 1) = 1;
                %   kzend = 64 - 1 = 63;
                % When kr = krmax / 2:  kzedge = (NIMM / lamda) * (1 - cosalpha) = 1.9298 um-1
                %   kzstart = round(1.9298 * 8 + 1) = 16
                %   kzend = 64 - 16 = 48
                for iz = kzstart : kzend
                    avg_output(ix + 1, iz + 1) = 0;
                end
                
            elseif order == 1
                jotfshape = round((kzedge / dkz) + 0.999);                
                kzend = ileave1_1 - jotfshape;          % ileave1_1 is the minimum shifted-pixel number for order1 OTF;
                % When kr = 0 or krmax: kzedge = (NIMM / lamda) * (cosalpha - cosalpha) = 0
                %   jotfshape = round(0 / dkz + 0.999) = 1;
                %   kzend = 16 - 1 = 15;
                % When kr = krmax / 2:  kzedge = (NIMM / lamda) * (1 - cosalpha) = 1.9298 um-1
                %   jotfshape = round(1.9298 * 8 + 0.999) = 16
                %   kzend = 16 - 16 = 0
                if kzend >= 0
                    for iz = 0 : kzend
                        if (iz == 0)
                            % iz + 1 = 1
                            avg_output(ix + 1, 1) = 0;
                        else
                            % iz + 1 = 2:16
                            avg_output(ix + 1, iz + 1) = 0;
                            % nz - iz + 1 = 50:64
                            avg_output(ix + 1, nz - iz + 1) = 0;
                        end
                    end
                end
                kzstart = ileave1_2 + jotfshape;
                % When kr = 0 or krmax: jotfshape = round(0 / dkz + 0.999) = 1
                %   kzstart = 20 + 1 = 21;
                % When kr = krmax / 2:  jotfshape = round(1.9298 * 8 + 0.999) = 16
                %   kzstart = 20 + 16 = 36
                for iz = kzstart : fix(nz / 2)
                    % When kr = 0 or krmax: iz + 1 = 22:33; nz - iz + 1 = 44:33
                    % That is 22:44
                    % When kr = krmax / 2:  (kzstart + 1) : fix(nz / 2) = 37:32, it is a empty double row vector
                    avg_output(ix + 1, iz + 1) = 0;
                    avg_output(ix + 1, nz - iz + 1) = 0;
                end               
                
            else
                % order == 2;
                jotfshape = round((kzedge / dkz) + 0.999);
                kzstart = ileave2 + jotfshape;          % ileave2 is the maximum axially-broadened pixel number for order2 OTF;
                % When kr = 0 or krmax: jotfshape = round(0 / dkz + 0.999) = 1
                %   kzstart = 3 + 1 = 4;
                % When kr = krmax / 2:  jotfshape = round(1.9298 * 8 + 0.999) = 16
                %   kzstart = 3 + 16 = 19;
                for iz = kzstart : fix(nz / 2)
                    % When kr = 0 or krmax: ix + 1 : nz - iz + 1 = 5 : 61
                    % When kr = krmax / 2:  ix + 1 : nz - iz + 1 = 20 : 46
                    avg_output(ix + 1, iz + 1) = 0;
                    avg_output(ix + 1, nz - iz + 1) = 0;
                end
            end
            
        % Clean up outside of lateral resolution limit
        % ix = 106:128
        else
            for iz = 0 : (nz - 1)
                avg_output(ix + 1, iz + 1) = 0;
            end
        end
        
    end
    
else
%     lambdaem = lamda / NIMM;
%     lambdaexc = 0.93 * lambdaem;                 % 0.93 approximates a typical lambdaexc/lambdaem
%     two_over_lambdaem = 2 / lambdaem;
%     two_over_lambdaexc = 2 / lambdaexc;
%     alpha = asin(NA_local / NIMM);              % aperture angle of objectives
%     beta = asin(k0mag / two_over_lambdaexc);    % angle of side illumination beams
%     betamin = asin((k0mag / two_over_lambdaexc) - sin(alpha) * SPOTRATIO);
%     
%     for ix = 0 : (icleanup - 1)
%         kr = ix * dkr;
%         if kr <= krmax
%             if order == 0
%                 center_of_arc = ceil(two_over_lambdaexc / dkz) + 3;
%             elseif order == 1
%                 center_of_arc = ceil((1 / lambdaexc) * (cos(beta) + 1) / dkz) + 3;
%             else
%                 center_of_arc = ceil(two_over_lambdaexc * cos(beta) / dkz) + 3;
%             end
%             kzstart = ceil(center_of_arc + sqrt(two_over_lambdaem * two_over_lambdaem - kr * kr) / dkz);
%             kzend = nz - kzstart;
%             for iz = kzstart : kzend
%                 avg_output(ix + 1, iz + 1) = 0;
%             end
%         else
%             for iz = 0 : (nz - 1)
%                 avg_output(ix + 1, iz + 1) = 0;
%             end
%         end
%     end

    for ix = 0 : (icleanup - 1)         % ix = 0:128
        kr = ix * dkr;                  % kr = 0:6.25
        % ix = 0:105
        if kr <= krmax
            beta = asin((NA_local - kr * lamda) / NIMM);
            % When kr = 0:          beta = asin(NA_local / NIMM) = asin(1.35 / 1.406) = 1.2876 (73.7748 degrees);
            % When kr = krmax / 2:  beta = 0
            % When kr = krmax:      beta = -1.2876 (-73.7748 degrees);
            % It can be regarded as: beta = asin((krmax / 2 - kr) / (NIMM / lamda));
            kzedge = (NIMM / lamda) * (cos(beta) - cosalpha);
%             fprintf("kzedge=%f, dkz=%f\n", kzedge, dkz);
            if order == 0
                kzstart = round((kzedge / dkz) + 1);
                jotfshape = round((kzedge / dkz) + 0.999);  
                kzend = ileave0_1 - jotfshape;
                if kzend >= kzstart
                    for iz = kzstart : kzend
                        if (iz == 0)
                            % iz + 1 = 1
                            avg_output(ix + 1, 1) = 0;
                        else
                            % iz + 1 = 2:16
                            avg_output(ix + 1, iz + 1) = 0;
                            % nz - iz + 1 = 50:64
                            avg_output(ix + 1, nz - iz + 1) = 0;
                        end
                    end
                end
                
                kzstart = ileave0_2 + jotfshape;
                % When kr = 0 or krmax: jotfshape = round(0 / dkz + 0.999) = 1
                %   kzstart = 20 + 1 = 21;
                % When kr = krmax / 2:  jotfshape = round(1.9298 * 8 + 0.999) = 16
                %   kzstart = 20 + 16 = 36
                for iz = kzstart : fix(nz / 2)
                    % When kr = 0 or krmax: iz + 1 = 22:33; nz - iz + 1 = 44:33
                    % That is 22:44
                    % When kr = krmax / 2:  (kzstart + 1) : fix(nz / 2) = 37:32, it is a empty double row vector
                    avg_output(ix + 1, iz + 1) = 0;
                    avg_output(ix + 1, nz - iz + 1) = 0;
                end  
                
            elseif order == 1
                jotfshape = round((kzedge / dkz) + 0.999);                
                kzend = ileave1_1 - jotfshape;          % ileave1_1 is the minimum shifted-pixel number for order1 OTF;
                % When kr = 0 or krmax: kzedge = (NIMM / lamda) * (cosalpha - cosalpha) = 0
                %   jotfshape = round(0 / dkz + 0.999) = 1;
                %   kzend = 16 - 1 = 15;
                % When kr = krmax / 2:  kzedge = (NIMM / lamda) * (1 - cosalpha) = 1.9298 um-1
                %   jotfshape = round(1.9298 * 8 + 0.999) = 16
                %   kzend = 16 - 16 = 0
                if kzend >= 0
                    for iz = 0 : kzend
                        if (iz == 0)
                            % iz + 1 = 1
                            avg_output(ix + 1, 1) = 0;
                        else
                            % iz + 1 = 2:16
                            avg_output(ix + 1, iz + 1) = 0;
                            % nz - iz + 1 = 50:64
                            avg_output(ix + 1, nz - iz + 1) = 0;
                        end
                    end
                end
                kzstart = ileave1_2 + jotfshape;
                jotfshape = round((kzedge / dkz) + 0.999);  
                kzend = ileave1_3 - jotfshape;
                % When kr = 0 or krmax: jotfshape = round(0 / dkz + 0.999) = 1
                %   kzstart = 20 + 1 = 21;
                % When kr = krmax / 2:  jotfshape = round(1.9298 * 8 + 0.999) = 16
                %   kzstart = 20 + 16 = 36
                if kzend >= kzstart
                    for iz = kzstart : kzend
                        if (iz == 0)
                            % iz + 1 = 1
                            avg_output(ix + 1, 1) = 0;
                        else
                            % iz + 1 = 2:16
                            avg_output(ix + 1, iz + 1) = 0;
                            % nz - iz + 1 = 50:64
                            avg_output(ix + 1, nz - iz + 1) = 0;
                        end
                    end
                end
                
                kzstart = ileave1_4 + jotfshape;
                % When kr = 0 or krmax: jotfshape = round(0 / dkz + 0.999) = 1
                %   kzstart = 20 + 1 = 21;
                % When kr = krmax / 2:  jotfshape = round(1.9298 * 8 + 0.999) = 16
                %   kzstart = 20 + 16 = 36
                for iz = kzstart : fix(nz / 2)
                    % When kr = 0 or krmax: iz + 1 = 22:33; nz - iz + 1 = 44:33
                    % That is 22:44
                    % When kr = krmax / 2:  (kzstart + 1) : fix(nz / 2) = 37:32, it is a empty double row vector
                    avg_output(ix + 1, iz + 1) = 0;
                    avg_output(ix + 1, nz - iz + 1) = 0;
                end               
                
            else
                % order == 2;
                jotfshape = round((kzedge / dkz) + 0.999);
                kzstart = ileave2 + jotfshape;          % ileave2 is the maximum axially-broadened pixel number for order2 OTF;
                % When kr = 0 or krmax: jotfshape = round(0 / dkz + 0.999) = 1
                %   kzstart = 3 + 1 = 4;
                % When kr = krmax / 2:  jotfshape = round(1.9298 * 8 + 0.999) = 16
                %   kzstart = 3 + 16 = 19;
                for iz = kzstart : fix(nz / 2)
                    % When kr = 0 or krmax: ix + 1 : nz - iz + 1 = 5 : 61
                    % When kr = krmax / 2:  ix + 1 : nz - iz + 1 = 20 : 46
                    avg_output(ix + 1, iz + 1) = 0;
                    avg_output(ix + 1, nz - iz + 1) = 0;
                end
            end
            
        % Clean up outside of lateral resolution limit
        % ix = 106:128
        else
            for iz = 0 : (nz - 1)
                avg_output(ix + 1, iz + 1) = 0;
            end
        end
        
    end
end
end

%% Fix the origin values using linear fitting
function avg_output = fixorigin(otfkxkz, nx, nz, kx1, kx2, dkr, dkz, lamdanm, twolens, NA, NIMM)
% Lin's notes:
% Theoretically the origin of a 3D OTF is a singularity
% -- it's value should be infinity!!
% Practically, a measured OTF will never have singularity at its origin,
% but it can be easily off by a lot especially because
% image background of all pixels all contribute to the OTF value at the origin.

% On the other hand, a 2D OTF's radial profile is almost linear near the origin
% and a 2D OTF is just a projection of 3D OTF along the kz axis.
% We therefore wanted to extrapolate 3D OTF's origin value from
% a line fit of its corresponding 2D OTF using
% a short segment of lateral pixels near, but not including, the origin.

% INPUT VARIABLES
%   otfkxkz = avg_output(:, :, phase);
%   nx = header.nx = 256;
%   nz = header.nz / nphases = 64;
%   kx1 = interpkr(1);          Suppose kx1 = 3;
%   kx2 = interpkr(2);          Suppose kx2 = 20;
%   kx1 < kx2;
% OUTPUT VARIABLE
%   avg_output

avg_output = otfkxkz;
totsum = 0;
ysum = 0;
sqsum = 0;

fprintf("In fixorigin()\n");
meani = 0.5 * (kx1 + kx2);
numvals = kx2 - kx1 + 1;
sum = zeros([kx2 + 1, 1], 'single');

% if (~twolens)
    for i = 0 : kx2
        sum(i + 1) = 0;
        
        if i == 0
            % don't want to add up garbages on kz axis
            sum(i + 1) = real(avg_output(1, 1));
        else
            for j = 0 : (nz - 1)
                % imaginary parts will cancel each other because of symmetry
                sum(i + 1) = sum(i + 1) + real(avg_output(i + 1, j + 1));
            end
        end
        
        if i >= kx1
            totsum = totsum + sum(i + 1);
            ysum = ysum + sum(i + 1) * (i - meani);
            sqsum = sqsum + (i - meani) * (i - meani);
        end
    end
    
% else
%     NA_local = NA;
%     lamda = lamdanm * 0.001;
%     sinalpha = NA_local / NIMM;
%     cosalpha = cos(asin(sinalpha));
%     
%     for i = 0 : kx2
%         sum(i + 1) = 0;
%         kr = i * dkr;
%         beta = asin((NA_local - kr * lamda) / NIMM);
%         kzedge = (NIMM / lamda) * (cos(beta) - cosalpha);
%         kzstart = round((kzedge / dkz) + 1);
%         kzend = nz - kzstart;
%         
%         if i == 0
%             % don't want to add up garbages on kz axis
%             sum(i + 1) = real(avg_output(1, 1));
%         else
%             for j = 0 : (kzstart - 1)
%                 % imaginary parts will cancel each other because of symmetry
%                 sum(i + 1) = sum(i + 1) + real(avg_output(i + 1, j + 1));
%             end
%             for j = (kzend + 1) : (nz - 1)
%                 % imaginary parts will cancel each other because of symmetry
%                 sum(i + 1) = sum(i + 1) + real(avg_output(i + 1, j + 1));
%             end
%         end
%         
%         if i >= kx1
%             totsum = totsum + sum(i + 1);
%             ysum = ysum + sum(i + 1) * (i - meani);
%             sqsum = sqsum + (i - meani) * (i - meani);
%         end
%     end
% end

slope = ysum / sqsum;
avg = totsum / numvals;

for i = 0 : (kx1 - 1)
    lineval = avg + (i - meani) * slope;
    avg_output(i + 1, 1) = avg_output(i + 1, 1) - (sum(i + 1) - lineval);
end

end

%% Scale the maximum value of order0 OTF to 1, and similar to order1 and order2
function [avg_output, scalefactor_return] = rescale(otfkxkz, order, nx, nz, scalefactor, dorescale)
% INPUT VARIABLES
%   otfkxkz = avg_output(:, :, phase);
%   order = fix(phase / 2) = 0, 1, 2;
%   nx = header.nx = 256;
%   nz = header.nz / nphases = 64;
%   scalefactor = 1.0
%   dorescale = 0 or 1;     Whether to do rescale for order0, order1 and order2 respectively
% OUTPUT VARIABLE
%   avg_output
%   scalefactor_return


% need to find out scalefactor
if (order == 0) || (dorescale)
    otf_mag = abs(otfkxkz);
    scalefactor = single(1 / max(otf_mag(:)));
end

avg_output = single(otfkxkz .* scalefactor);
scalefactor_return = scalefactor;

end

%% Find out the phase of the side bands relative to the center band and get the real OTF of the side bands.
function otf_return = combine_reim(otf, norders, nx, nz, bForcedPIshift)
% It's based on the fact that
%   bandre = OTF * cos(phi), and bandim = OTF * sin(-phi)
% So we have: phi = atan(abs(bandim) / abs(bandre));
% In most cases, OTF is a complex number with real and imag part are <both positive>.
% Therefor, there are 4 cases:
%   1. real(bandre) > 0 and real(bandim) > 0:   phi = -phi;
%   2. real(bandre) > 0 and real(bandim) < 0:   phi = phi;
%   3. real(bandre) < 0 and real(bandim) > 0:   phi = phi + pi;
%   4. real(bandre) < 0 and real(bandim) < 0:   phi = pi - phi;
% With this definition, bandplus = bandre + i bandim
%                       bandminus = bandre - i bandim
% After we calculated phi, we have:
%   OTF = bandre * cos(phi) + bandim * sin(-phi) = OTF * cos(phi) ^ 2 + OTF * sin(phi) ^ 2 = OTF

% It should be noted that, the "phi" here is inverse to which defined in Mats' 2008 Biophysical Journal paper
% In his paper:
%   bandre = OTF * cos(phi), and bandim = OTF * sin(phi)
%   bandplus = bandre + i bandim
%   bandminus = bandre - i bandim
% It’s just a matter of signs

% INPUT VARIABLES
%   otf = avg_output(:, :, :);
%   norders = (nphases + 1) / 2 = 3;
%   nx = header.nx = 256;
%   nz = header.nz / nphases = 64;
%   bForcedPIshift;                     Whether to force PI change
% OUTPUT VARIABLE
%   avg_output
%   scalefactor_return

otf_return = otf;

fprintf("In combine_reim()\n");
index_x = fix(nx / 8);          % index_x = 32


for order = 1 : (norders - 1)
    bandre_mag = sum(abs(otf(:, :, 2 * order)), 'all');
    bandim_mag = sum(abs(otf(:, :, 2 * order + 1)), 'all');
    
    phi = atan(bandim_mag / bandre_mag);
    % phi is now in the first quadrant only, since bandim_mag and bandre_mag are both positive
    % which quadrant phi should be in is decided by the kz = 0 plane of bandre and bandim values

    if (real(otf(index_x + 1, 1, 2 * order)) < 0) && (real(otf(index_x + 1, 1, 2 * order + 1)) > 0)
        % 3rd quadrant
        phi = phi + pi;
    elseif (real(otf(index_x + 1, 1, 2 * order)) < 0) && (real(otf(index_x + 1, 1, 2 * order + 1)) < 0)
        % 2nd quadrant
        phi = pi - phi;
    elseif (real(otf(index_x + 1, 1, 2 * order)) > 0) && (real(otf(index_x + 1, 1, 2 * order + 1)) > 0)
        % 4th quadrant
        phi = -phi;
    end
    
    % Sometimes the above logic still won't find the correct phase; user can specify additional Pi phase shift
    if (order == 1) && (bForcedPIshift)
        phi = pi - phi;
    end
    
    fprintf("  phi=%f\n", phi);
    
    otf_return(:, :, 2 * order) = otf(:, :, 2 * order) .* cos(phi) + otf(:, :, 2 * order + 1) .* sin(-phi);
    otf_return(:, :, 2 * order + 1) = 0;   % since bandim = band * sin(phi) and phi is essentially 0;

end
end
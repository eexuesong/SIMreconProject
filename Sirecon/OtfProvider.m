classdef OtfProvider
    properties
        otffiles = [];
        norders = [];
        bRadAvgOTF = [];
        bOneOTFperAngle = [];
        nDirsOTF = [];
        bTwolens = [];
        
        nxotf = [];     % OTF's x pixel number
        nyotf = [];     % OTF's y pixel number
        nzotf = [];     % OTF's z pixel number
        dkrotf = [];    % Radial-direction inverse microns per pixel in OTF
        dkzotf = [];    % Axial-direction inverse microns per pixel in OTF
        
        vals = [];
        valsAtt = [];
        valsOnlyAtt = [];
        
        na = [];
        lamda = [];
        krcutoff = [];
        rdistcutoff = [];
        
    end
    
    methods (Access = public)
        function obj = OtfProvider(params, imgParams)
            
            obj.otffiles = params.otffiles;
            obj.norders = params.norders;
            obj.bRadAvgOTF = params.bRadAvgOTF;
            obj.bOneOTFperAngle = params.bOneOTFperAngle;
            if (obj.bOneOTFperAngle)
                obj.nDirsOTF = params.ndirs;
            else
                obj.nDirsOTF = 1;
            end
            obj.bTwolens = params.bTwolens;
            
            obj.na = params.na;
            obj.lamda = imgParams.wave(1);
            % "krcutoff" is the lateral maximum frequency of emmision wavelength in current OTF.
            obj.krcutoff = 2.0 * obj.na / (obj.lamda / 1.0e3); 
            
        end
        
        %% To get the OTF header information from otf-file. So (nxotf, nyotf, nzotf, dkrotf, dkzotf) are stored in m_myParams.
        function obj = determine_otf_dimensions(obj, imgParams)
            
            [otfstream_no, ~] = fopen(obj.otffiles, 'r');
            otfheader = ImageJ_formatted_TIFF.parse_tif(otfstream_no, 0);
            
            % determine nzotf, nxotf, nyotf, dkrotf, dkzotf based on dataset being 2D/3D and flags bRadAvgOTF and bOneOTFperAngle
            if imgParams.nz == 1
                % 2D, ignore bOneOTFperAngle
                obj.nxotf = otfheader.ImageWidth;
                if (obj.bRadAvgOTF)
                    obj.nyotf = 1;
                    obj.nxotf = 2 * (otfheader.ImageLength - 1);        % Avoid decreasing x by 2 due to radial average, default ImageWidth is a even.
                else
                    obj.nxotf = otfheader.ImageWidth;
                    obj.nyotf = otfheader.ImageLength;
                end
                obj.nzotf = 1;
                obj.dkrotf = 1 / (otfheader.resolution * obj.nxotf);    % dkrotf's unit is 1/micron
                obj.dkzotf = 1;
            else
                % 3D
                % otfheader.nz is defined as the total number of sections = nz * (ndirs * nphases) in our case
                ImageDescription = convertStringsToChars(otfheader.ImageDescription);
                % Calculate nz from header.ImageDescription
                k1 = strfind(ImageDescription, 'images=');
                ImageDescription_crop = ImageDescription(k1:end);
                k2 = strfind(ImageDescription_crop, newline);
                if ~isempty(k1) && ~isempty(k2)
                    otfheadernz = str2double(ImageDescription(k1(1) + 7: k1(1) + k2(1) - 2));
                else
                    warning("Did not find 'images=' in ImageDescription. Try to calculate it from filesize.");
                    FileAttributes = dir(otffilein);
                    otfheadernz = floor(FileAttributes.bytes / (otfheader.ImageWidth * otfheader.ImageLength * (otfheader.BitsPerSample / 8)));
                end
                
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
                
                if (obj.bRadAvgOTF)
                    obj.nzotf = otfheader.ImageWidth;
                    obj.nxotf = 2 * (otfheader.ImageLength - 1);      % Avoid decreasing x by 2 due to radial average, default ImageWidth is a even.
                    obj.nyotf = 1;
                    obj.dkzotf = 1 / (spacing * obj.nzotf);               % dkzotf's unit is 1/micron
                    obj.dkrotf = 1 / (otfheader.resolution * obj.nxotf);  % dkrotf's unit is 1/micron
                else
                    % each order has a 3D OTF stack (non-negative kx half of Fourier space)
                    obj.nzotf = fix(otfheadernz / obj.norders);
                    if (obj.bOneOTFperAngle)
                        obj.nzotf = fix(obj.nzotf / obj.ndirs);
                    end
                    obj.nxotf = otfheader.ImageWidth;
                    obj.nyotf = otfheader.ImageLength;
                    obj.dkzotf = 1 / (spacing * obj.nzotf);
                    obj.dkrotf = 1 / (otfheader.resolution * obj.nxotf);
                    
                end
            end
            % "rdistcutoff" is the lateral pixels of maximum frequency of emmision wavelength in current OTF.
            obj.rdistcutoff = obj.krcutoff / obj.dkrotf;
            fclose(otfstream_no);
            fprintf("nzotf=%d, dkzotf=%f um-1, nxotf=%d, nyotf=%d, dkrotf=%fum-1\n", obj.nzotf, obj.dkzotf, obj.nxotf, obj.nyotf, obj.dkrotf);
        end        
        
        %% Transfer OTF data from file to RAM
        function obj = loadOTFs(obj)
            
            otfmat = strrep(obj.otffiles, 'tif', 'mat');
            load(otfmat, 'avg_output_mat');
            
            % If radially-averaged OTF is used, obj.nyotf = 1
            % size(avg_output_mat) = [fix(obj.nxotf / 2) + 1, obj.nzotf, obj.norders, obj.nDirsOTF]
            if (obj.bRadAvgOTF)
                % padding "rxy" dimension from fix(obj.nxotf / 2) + 1 to obj.nxotf
                avg_output_mat(end + 1 : obj.nxotf, :, :, :) =  complex(single(0.0), single(0.0));
                if (obj.nzotf == 1)
                    % 2D OTF is just a projection of 3D OTF along the kz axis.
                    avg_output_mat = sum(avg_output_mat, 2);
                end
            else
                if (obj.nzotf == 1)
                    % 2D OTF is just a projection of 3D OTF along the kz axis.
                    avg_output_mat = sum(avg_output_mat, 3);
                end
            end
            
            obj.vals = complex(zeros([obj.nyotf, obj.nxotf, obj.nzotf, obj.norders, obj.nDirsOTF], 'single'), zeros([obj.nyotf, obj.nxotf, obj.nzotf, obj.norders, obj.nDirsOTF], 'single'));
            try
                obj.vals(:) = avg_output_mat(:);
            catch
                error('OTF .tif and .mat files have different dimensions!\n');
            end
        end
        
        %% Return the OTF cutoff, unit is 
        function cutoff = getCutoff(obj)
            cutoff = obj.krcutoff;
        end
        
        %% Create a 3-dimension, radial symmetric matrix of the OTF, centered at kx,ky. Attenuation is not applied. 
        function otfval = getOtfMatrix3D(obj, order, dir, nx, ny, kx, ky, krscale, nz, kzscale)
            %{
                order:  OTF order
                dir:    OTF direction
                nx:     X pixel number of output matrix
                ny:     Y pixel number of output maxtix
                kx:     Position / offset kx
                ky:     Poistion / offset ky
                krscale:the ratio between OTF FOV and image FOV, latteraly;
                nz:     Z pixel number of output matrix, (2 * nz + 1) due to symmetry
                kzscale:the ratio between OTF FOV and image FOV, axially;
            %}
            if ((order < 0) || (order >= obj.norders))
                error("Order index exceeds OTF imported or <0")
            end
            
            if ((dir < 1) || (dir > obj.nDirsOTF))
                error("Direction index exceeds OTF imported or <1")
            end
            
            if obj.bRadAvgOTF
                x = 0 : (nx - 1);
                y = 0 : (ny - 1);
                z = -nz : nz;
                x(x >= (nx / 2)) = x(x >= (nx / 2)) - nx;
                y(y >= (ny / 2)) = y(y >= (ny / 2)) - ny;
                x = x + kx;
                y = y + ky;
                [X, Y, Z] = meshgrid(x, y, z);
                % Because OTF is an even function, we only need to know the lateral distance from each pixel to the center point, and neglect (x, y);
                krindex = sqrt(X .* X + Y .* Y) .* krscale;
                % Avoid one situation that: krindex value is larger than fix(obj.nxotf / 2)
                krindex(krindex > fix(obj.nxotf / 2)) = fix(obj.nxotf / 2) - 1;
                
                % kzindex now is the axial distance from each pixel to the center point 0 in OTF after fftshift;
                kzindex = Z .* kzscale;
                % kzindex now is the axial pixel to the center point 0 in OTF before fftshift;
                % Suppose nzotf = 100, then z = [-20:20] --> [80:99]...[0:20] = kzindex;
                kzindex(kzindex < 0) = kzindex(kzindex < 0) + obj.nzotf;
                
                irindex = floor(krindex);
                izindex = floor(kzindex);
                
                ar = krindex - irindex;     % decimal part of lateral index
                az = kzindex - izindex;     % decimal part of axial index
                % "Four-corner" indices: normal situation
                upper_left = irindex .* obj.nzotf + izindex;
                upper_right = irindex .* obj.nzotf + (izindex + 1);
                lower_left = (irindex + 1) .* obj.nzotf + izindex;
                lower_right = (irindex + 1) .* obj.nzotf + (izindex + 1);
                
                % Boundary conditions: right side
                upper_right(izindex == obj.nzotf - 1) = irindex(izindex == obj.nzotf - 1) .* obj.nzotf;
                lower_right(izindex == obj.nzotf - 1) = (irindex(izindex == obj.nzotf - 1) + 1) .* obj.nzotf ;
                
                % Change otf from [1, r, z] order into [z, 1, r] order
                otf_swap = permute(obj.vals(:, :, :, order + 1, dir), [3, 1, 2]);
                column_otf = otf_swap(:);
                ar = ar(:);
                az = az(:);
                upper_left = upper_left(:);
                upper_right = upper_right(:);
                lower_left = lower_left(:);
                lower_right = lower_right(:);
                
                % otfval now is a column-wize vector
                otfval = (1 - ar) .* (column_otf(upper_left + 1) .* (1 - az) + column_otf(upper_right + 1) .* az) + ...
                    ar .* (column_otf(lower_left + 1) .* (1 - az) + column_otf(lower_right + 1) .* az);
                % reshape otfval back to its correct size
                otfval = reshape(otfval, [ny, nx, (2 * nz + 1)]);
                % clear the region where krindex value is out of OTF support (obj.rdistcutoff)
                otfval(krindex > obj.rdistcutoff) = 0; 
                % clear the region where krindex value is out of matrix support fix(obj.nxotf / 2)
                otfval(krindex > fix(obj.nxotf / 2)) = 0;
            else
                error("Sorry, this program only handles radial-averaged OTF\n")
            end
            
        end
        
        %% Create a 2-dimension, radial symmetric matrix of the OTF, centered at kx,ky. Attenuation is not applied. 
        function otfval = getOtfMatrix2D(obj, order, dir, nx, ny, kx, ky, krscale)
            %{
                order:  OTF order
                dir:    OTF direction
                nx:     X pixel number of output matrix
                ny:     Y pixel number of output maxtix
                kx:     Position / offset kx
                ky:     Poistion / offset ky
                krscale:the ratio between OTF FOV and image FOV, latteraly;
            %}
            if ((order < 0) || (order >= obj.norders))
                error("Order index exceeds OTF imported or <0")
            end
            
            if ((dir < 1) || (dir > obj.nDirsOTF))
                error("Direction index exceeds OTF imported or <1")
            end
            
            if obj.bRadAvgOTF
                x = 0 : (nx - 1);
                y = 0 : (ny - 1);
                x(x >= (nx / 2)) = x(x >= (nx / 2)) - nx;
                y(y >= (ny / 2)) = y(y >= (ny / 2)) - ny;
                x = x + kx;
                y = y + ky;
                [X, Y] = meshgrid(x, y);
                % Because OTF is an even function, we only need to know the lateral distance from each pixel to the center point, and neglect (x, y);
                krindex = sqrt(X .* X + Y .* Y) .* krscale;
                % Avoid one situation that: krindex value is larger than fix(obj.nxotf / 2)
                krindex(krindex > fix(obj.nxotf / 2)) = fix(obj.nxotf / 2) - 1;
                irindex = floor(krindex);
                ar = krindex - irindex;     % decimal part of lateral index
                
                upper_index = irindex;
                lower_index = irindex + 1;
                
                if (obj.nzotf == 1)
                    % Already a 2D OTF
                    otf = obj.vals(:, :, :, order + 1, dir);
                else
                    % 2D OTF is just a projection of 3D OTF along the kz axis.
                    otf = sum(obj.vals(:, :, :, order + 1, dir), 3);
                end
                
                % Change otf from [1, r, z] order into [z, 1, r] order
                otf_swap = permute(otf, [3, 1, 2]);
                column_otf = otf_swap(:);
                ar = ar(:);
                upper_index = upper_index(:);
                lower_index = lower_index(:);
                
                % otfval now is a column-wize vector
                otfval = (1 - ar) .* (column_otf(upper_index + 1)) + ar .* (column_otf(lower_index + 1));
                % reshape otfval back to its correct size
                otfval = reshape(otfval, [ny, nx]);
                % clear the region where krindex value is out of OTF support (obj.rdistcutoff)
                otfval(krindex > obj.rdistcutoff) = 0; 
                % clear the region where krindex value is out of matrix support fix(obj.nxotf / 2)
                otfval(krindex > fix(obj.nxotf / 2)) = 0;
            else
                error("Sorry, this program only handles radial-averaged OTF\n")
            end
            
        end
        
        %% Creates an attenuation matrix (for optical sectioning).
        function attval = getAttenuationMatrix2D(obj, strength, fwhm, nx, ny, kx, ky, krscale)
            x = 0 : (nx - 1);
            y = 0 : (ny - 1);
            x(x >= (nx / 2)) = x(x >= (nx / 2)) - nx;
            y(y >= (ny / 2)) = y(y >= (ny / 2)) - ny;
            x = x + kx;
            y = y + ky;
            [X, Y] = meshgrid(x, y);            
            krindex = sqrt(X .* X + Y .* Y) .* krscale;
            
            cycl = krindex ./ obj.dkrotf;
            attval = obj.valAttenuation(cycl, strength, fwhm);
        end
        
        %% Creates an 3D attenuation matrix.
        function attval = getAttenuationMatrix3D(obj, order, nx, ny, kx, ky, nz, rdistcutoff, params)
            
            if ((order < 0) || (order >= obj.norders))
                error("Order index exceeds OTF imported or <0")
            end
            
            x = 0 : (nx - 1);
            y = 0 : (ny - 1);
            z = -nz : nz;
            x(x >= (nx / 2)) = x(x >= (nx / 2)) - nx;
            y(y >= (ny / 2)) = y(y >= (ny / 2)) - ny;
            x = x + kx;
            y = y + ky;
            [X, Y, Z] = meshgrid(x, y, z);
            rdist = sqrt(X .* X + Y .* Y);
            
            if (params.bSuppress_singularities && order ~= 0)
                attval = obj.suppress_singularities(rdist, params.dampenFactor);
                attval(rdist > params.suppression_radius) = 1;
            elseif (~params.bDampenOrder0 && params.bSuppress_singularities && order == 0)
                attval = obj.suppress_singularities(rdist, params.dampenFactor);
                attval(rdist > params.suppression_radius) = 1;
            elseif (params.bDampenOrder0 && order == 0)
                attval = obj.order0damping(rdist, Z, rdistcutoff, nz);
            end
            
            % if no kz=0 plane is used:
            if (order == 0 && params.bNoKz0)
                attval(Z == 0) = 0;
            end
            
            attval(rdist > rdistcutoff) = 0;
        end
        
        %% Set components not covered by OTF support to zero.
        function mat = maskOtf(obj, mat, kx, ky, rdistcutoff)
            
            [ny, nx, nz] = size(mat);
            x = 0 : (nx - 1);
            y = 0 : (ny - 1);
            z = 0 : (nz - 1);
            x(x >= (nx / 2)) = x(x >= (nx / 2)) - nx;
            y(y >= (ny / 2)) = y(y >= (ny / 2)) - ny;
            x = x + kx;
            y = y + ky;
            [X, Y, ~] = meshgrid(x, y, z);
            rdist = sqrt(X .* X + Y .* Y);
            mat(rdist > rdistcutoff) = 0;
            
        end
        
        %% Create an apodization matrix.
        function apofact = writeApoFactor(obj, nx, ny, kx, ky, nz, apocutoff, zapocutoff, params)
            x = 0 : (nx - 1);
            y = 0 : (ny - 1);
            z = -nz : nz;
            x(x >= (nx / 2)) = x(x >= (nx / 2)) - nx;
            y(y >= (ny / 2)) = y(y >= (ny / 2)) - ny;
            x = x + kx;
            y = y + ky;
            [X, Y, Z] = meshgrid(x, y, z);
            rdist = sqrt(X .* X + Y .* Y);
            
            zdistabs = abs(Z);
            
            if (zapocutoff > 0)
                % 3D case
                rho = sqrt((rdist ./ apocutoff) .^ 2 + (zdistabs ./ zapocutoff) .^ 2);
            else
                % 2D case
                rho = sqrt((rdist ./ apocutoff) .^ 2);
            end
            
            rho(rho > 1.0) = 1.0;
            
            if (params.apodizeoutput == 1)
                % cosine-apodize
                apofact = cos(0.5 * pi * rho);
            elseif (params.apodizeoutput == 2)
                % triangle-apodize
                apofact = 1.0 - rho .^ params.apoGamma;                
            end
        end
        
        %% Create a ButterWorth apodization matrix.
        function Butterworth_apofact = writeApoFactor_Butterworth(obj, nx, ny, kx, ky, nz, apocutoff, params)
            % apocutoff : width of Butterworth Filter
            x = 0 : (nx - 1);
            y = 0 : (ny - 1);
            z = -nz : nz;
            x(x >= (nx / 2)) = x(x >= (nx / 2)) - nx;
            y(y >= (ny / 2)) = y(y >= (ny / 2)) - ny;
            x = x + kx;
            y = y + ky;
            [X, Y, ~] = meshgrid(x, y, z);
            rdist = sqrt(X .* X + Y .* Y);
            
            beta = 0.0001;  % Value at cut off pixel
            ee = 1 / beta^2 - 1;
            n = 50;         % Slope 
            Butterworth_apofact = 1 / sqrt(1 + ee .* ((rdist /  apocutoff).^2 .^ n));
            
        end
        
        %% Create a Gaussian apodization matrix.
        function Gaussian_apofact = writeGaussianApoFactor(obj, nx, ny, kx, ky, nz, apocutoff, zapocutoff, params)
            x = 0 : (nx - 1);
            y = 0 : (ny - 1);
            z = -nz : nz;
            x(x >= (nx / 2)) = x(x >= (nx / 2)) - nx;
            y(y >= (ny / 2)) = y(y >= (ny / 2)) - ny;
            x = x + kx;
            y = y + ky;
            [X, Y, Z] = meshgrid(x, y, z);
            rdist = sqrt(X .* X + Y .* Y);
            
            zdistabs = abs(Z);
            sigma_xy = apocutoff / sqrt(2 * log(200));
            sigma_z = zapocutoff / sqrt(2 * log(200));
            
            if (zapocutoff > 0)
                % 3D case
                Gaussian_apofact = exp(-rdist .^ 2 / (2 * sigma_xy ^ 2) - zdistabs .^ 2 / (2 * sigma_z ^ 2));
            else
                % 2D case
                Gaussian_apofact = exp(-rdist .^ 2 / (2 * sigma_xy ^ 2));
            end
            
%             Gaussian_apofact(rdist > apocutoff) = 0;
        end
        
        %% Create a 1D-Gaussian apodization matrix.
        function Butterworth_apofact = writeApoFactor_Butterworth_1D(obj, nx, ny, kx, ky, nz, apocutoff, params, theta)
            % apocutoff : width of Butterworth Filter
            x = 0 : (nx - 1);
            y = 0 : (ny - 1);
            z = -nz : nz;
            x(x >= (nx / 2)) = x(x >= (nx / 2)) - nx;
            y(y >= (ny / 2)) = y(y >= (ny / 2)) - ny;
            x = x + kx;
            y = y + ky;
            [X, Y, ~] = meshgrid(x, y, z);
            XY = [X(:) Y(:)];
            R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            rotXY = XY * R';
            Xr = reshape(rotXY(:,1), size(X,1), size(X,2), []);
            Yr = reshape(rotXY(:,2), size(Y,1), size(Y,2), []);
            
            beta = 0.0001;  % Value at cut off pixel
            ee = 1 / beta^2 - 1;
            n = 50;         % Slope 
            Butterworth_apofact = 1 / sqrt(1 + ee .* ((Yr /  apocutoff).^2 .^ n));
        end       

        
    end
        
    methods (Access = private)
        %% to dampen the OTF values near band centers using one minus a Gaussian of maximal amplitude (0...1) and fwhm
        function out = valAttenuation(obj, dist, strength, fwhm)
            % dist:     Distance to center in cycles/micron
            % strength: Strength of the attenuation, typically 0.95...0.999
            % fwhm:     FWHM of Attenuation, in cycles/micron, typically 0.8...2.0
            
            % a Gaussian function is in form:
            % f(x) = a * exp(-(x - b) ^ 2 / (2 * c ^ 2));
            % The parameter c is related to the full width at half maximum (FWHM) of the peak according to
            % FWHM = 2 * sqrt(2 * ln(2)) * c = 2.35482 * c;
            % So in our case, c = fwhm / 2.355
            
            % For NA = 1.35, lamda = 0.525 um, the maxmimum of OTF support is
            % 2 * NA / lamda  = 2 * 1.35 / 0.525 = 5.143 um -1
            % In real application, we do not want fwhm to be that high
            
            out = 1 - strength * exp((- dist .^ 2) / (2 * (fwhm / 2.355) ^ 2));
            
        end
        
        %% to dampen the OTF values near band centers;
        function out = suppress_singularities(obj, x, dampenFactor)
            
            x6 = x .^ 6;
            out = 1.0 ./ (1 + (20 * dampenFactor) ./ (x6 + 20));
            %{
                if x1 = 0; y1 = 0; rdist1 = 0; out = 1 / 1001 = 9.9900e-04;
                if x1 = 6; y1 = 8; rdist1 = 10; out = 0.9804;
                if x1 = 100; y1 = 100; rdist1 = 141.42; out = 1;
            %}
            
        end
        
        %% to dampen the OTF values near band centers for order = 0;
        function out =  order0damping(obj, radius, zindex, rlimit, zlimit)
            % radius = rdist1
            % zindex = Z0
            % rlimit = rdistcutoff
            % zlimit = zdistcutoff(1)
            
            rfraction = radius ./ rlimit;
            if (zlimit > 0)
                % 3D case
                zfraction = abs(zindex ./ zlimit);
            else
                % 2D case
                zfraction = 0;
            end
            
            out = rfraction .^ 2 + zfraction .^ 3;
            out(out > 1) = 1;
            
        end
    end
end
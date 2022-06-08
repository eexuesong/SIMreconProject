classdef ImageJ_formatted_TIFF
    methods(Static)
        function WriteTifStack(Stack, Filename)
            fileID = fopen(Filename, 'w+');
            [Height, Width, Depth] = size(Stack);
            minval = min(Stack, [], 'all');
            maxval = max(Stack, [], 'all');
            header = ImageJ_formatted_TIFF.write_IFD(fileID, Width, Height, class(Stack), 1, Depth, 1, 0.04, minval, maxval, 0.039);
            Stack = permute(Stack, [2 1 3]);
            Stack = reshape(Stack, [1, Height * Width * Depth]);
            
            [~, ~, system_endian] = computer;
            if system_endian ~= header.endian
                Stack = swapbytes(Stack);
            end
            Stack = typecast(Stack, 'uint8');
            fwrite(fileID, Stack, 'uint8');        
            fclose(fileID);
        end
        
        function Stack = ReadTifStack(Filename)
            % Read and parse tiff header first
            fileID = fopen(Filename, 'r');
            header = ImageJ_formatted_TIFF.parse_tif(fileID, 0);
            ImageWidth = header.ImageWidth;
            ImageLength = header.ImageLength;           
            ImageDescription = convertStringsToChars(header.ImageDescription);
            
            % Calculate slices from header.ImageDescription
            k1 = strfind(ImageDescription, 'images=');
            ImageDescription_crop = ImageDescription(k1:end);
            k2 = strfind(ImageDescription_crop, newline);
            if ~isempty(k1) && ~isempty(k2)
                Depth = str2double(ImageDescription(k1(1) + 7: k1(1) + k2(1) - 2));
            else
                warning("Did not find 'images=' in ImageDescription. Try to calculate it from filesize.");
                FileAttributes = dir(Filename);
                Depth = floor(FileAttributes.bytes / (ImageWidth * ImageLength * (header.BitsPerSample / 8)));
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
            PixelNum = ImageWidth * ImageLength;
            ByteCounts = fix(PixelNum * header.BitsPerSample / 8);
            Stack = zeros([1, Depth * PixelNum], dtype);
            
            fseek(fileID, header.StripOffsets, 'bof');
            for i = 1:Depth
                Stack((i - 1) * PixelNum + 1 : i * PixelNum) = typecast(transpose(fread(fileID, ByteCounts, 'uint8=>uint8')), dtype);
            end
            [~, ~, system_endian] = computer;
            if system_endian ~= header.endian
                Stack = swapbytes(Stack);
            end

            Stack = reshape(Stack, [PixelNum, 1, Depth]);
            Stack = reshape(Stack, [ImageWidth, ImageLength, Depth]);
            Stack = permute(Stack, [2 1 3]); 
            fclose(fileID);
        end
        
        function ifd = write_IFD(fileID, width, height, data_type, channels, slices, frames, spacing, min, max, resolution)
            %         We'll structure our TIF in the following way (same as how ImageJ saves tiff stack):
            %         8-byte Image File Header
            %         1st image file directory (IFD)
            %         Image description (~100 bytes)
            %         Image XResolution (Two 32-bit unsigned integers, 8 bytes)
            %         Image YResolution (Two 32-bit unsigned integers, 8 bytes)
            %         1st image data
            %         2nd image data
            %         ...
            %         last image data
            
            endian = 'L';   % I suggest to use endian = 'L'; because big-endian is slower because of swapbytes.
            [~, ~, system_endian] = computer;
            if system_endian == endian
                swapbytes_flag = 0;
            else
                swapbytes_flag = 1;
            end
            
            %% Write Tiff Header into file
            if endian == 'L'
                fwrite(fileID, sprintf("\x49\x49\x2A\x00"), 'uint8');   % little-endian (Intel format) order
            else
                fwrite(fileID, sprintf("\x4D\x4D\x00\x2A"), 'uint8');   % big-endian (Motorola format) order
            end
            IFDOffset_uint32 = 8;
            if swapbytes_flag == 0
                IFDOffset = typecast(uint32(IFDOffset_uint32), 'uint8');
            else
                IFDOffset = typecast(swapbytes(uint32(IFDOffset_uint32)), 'uint8');
            end
            fwrite(fileID, IFDOffset, 'uint8');
            
            %% IFD common part
            ifd = Simple_IFD(endian);
            ifd.ImageWidth = width;
            ifd.ImageLength = height;
            ifd.RowsPerStrip = height;
            ifd = ifd.set_dtype(data_type);
            ifd.StripByteCounts = fix(width * height * ifd.BitsPerSample / 8);
            ifd.NextIFD = 0;
            
            %% Image description part
            hyperstack_flag = 0;
            images = channels * slices * frames;
            image_description =  sprintf("ImageJ=1.53f\nimages=%d\n", images);
            if channels > 1
                hyperstack_flag = hyperstack_flag + 1;
                image_description = image_description + sprintf("channels=%d\n", channels);
            end            
            if slices > 1
                hyperstack_flag = hyperstack_flag + 1;
                image_description = image_description + sprintf("slices=%d\n", slices);
            end            
            if frames > 1
                hyperstack_flag = hyperstack_flag + 1;
                image_description = image_description + sprintf("frames=%d\n", frames);
            end            
            if hyperstack_flag > 1
                image_description = image_description + sprintf("hyperstack=true\n");
            end
            image_description = image_description + sprintf("mode=grayscale\n");
            % "\u00B5" is Unicode Character 'MICRO SIGN' (U+00B5)
            % Character '\u' is not valid in Matlab, so we use '\\u00B5' or '\x00B5'
            image_description = image_description + sprintf("unit=\\u00B5m\n");
            image_description = image_description + sprintf("spacing=%0.3f\n", spacing);
            image_description = image_description + sprintf("loop=false\n");
            image_description = image_description + sprintf("min=%0.1f\n", min);
            image_description = image_description + sprintf("max=%0.1f\n", max);
            image_description = image_description + sprintf("\x00");
            ifd.ImageDescription = image_description;
            ifd.NumberCharsInImageDescription = strlength(image_description);
            ifd.OffsetOfImageDescription = IFDOffset_uint32 + size(ifd.bytes, 2);
            
            %% Image resolution part
            image_resolution = zeros([1, 8], 'uint8');
            resolution_numerator = 1000000;       % Convert mm into nm. Although we use ifd.ResolutionUnit = 3 (Centimeter), ImageJ read unit from Metadata instead of ifd.ResolutionUnit.
            resolution_denominator = uint32(1000000 * resolution);    % nm / pixel, e.g. 39 nm
            
            if swapbytes_flag == 0
                image_resolution(1:4) = typecast(uint32(resolution_numerator), 'uint8');
                image_resolution(5:8) = typecast(uint32(resolution_denominator), 'uint8');
            else
                image_resolution(1:4) = typecast(swapbytes(uint32(resolution_numerator)), 'uint8');
                image_resolution(5:8) = typecast(swapbytes(uint32(resolution_denominator)), 'uint8');
            end
            ifd.XResolution = ifd.OffsetOfImageDescription + ifd.NumberCharsInImageDescription;
            ifd.YResolution = ifd.XResolution + 8;
            ifd.resolution =  round(single(resolution_denominator) / single(resolution_numerator), 5); % Unit: um / pixel
            
            %% Write IFD, ImageDescription and XYResolution into file
            ifd.StripOffsets = ifd.YResolution + 8;
            fwrite(fileID, ifd.bytes, 'uint8');
            fwrite(fileID, image_description, 'uint8');
            fwrite(fileID, image_resolution, 'uint8');
            fwrite(fileID, image_resolution, 'uint8');
        end
        
        function header = parse_tif(fileID, verbose)
            %%
            %             Open a file, determine that it's a TIF by parsing its header, and
            %             read through the TIF's Image File Directories (IFDs) one at a time
            %             to determine the structure of the TIF.
            %             See:
            %             partners.adobe.com/public/developer/en/tiff/TIFF6.pdf
            %             for reference.
            if nargin == 1
                verbose = 0;
            elseif nargin == 2
                verbose = verbose;
            end
            header = Simple_IFD();
            [next_ifd_offset, header.endian] = ImageJ_formatted_TIFF.parse_header(fileID, verbose);
            header = ImageJ_formatted_TIFF.parse_ifd(fileID, header, next_ifd_offset, verbose);
            
            
        end
        
        function [IFDOffset, endian] = parse_header(fileID, verbose)
            header = ImageJ_formatted_TIFF.get_bytes_from_file(fileID, 0, 8);
            if verbose
                fprintf("Header: ");
                for i = 1:size(header, 1)
                    fprintf("%d, ", header(i));
                end
                fprintf("\n");
            end
            if header(1:4) == [73; 73; 42; 0]
                endian = "L";
            elseif header(1:4) == [77; 77; 0; 42]
                endian = "B";
            else
                ME = MException("Not a TIF file");
                throw(ME)
            end
            IFDOffset = ImageJ_formatted_TIFF.bytes_to_int(header(5:8), endian);
            if verbose
                fprintf("This file is a %s-endian tif.\n", endian);
                fprintf("The offset of the first IFD is at %d bytes.\n", IFDOffset);
            end
        end
        
        function header = parse_ifd(fileID, header, ifd_offset, verbose)
            num_tags = ImageJ_formatted_TIFF.bytes_to_int(ImageJ_formatted_TIFF.get_bytes_from_file(fileID, ifd_offset, 2), header.endian);
            if verbose
                fprintf("IFD at offset %d bytes with %d tags:\n", ifd_offset, num_tags);
            end
            header.NumDirEntries = num_tags;
            ifd_bytes = ImageJ_formatted_TIFF.get_bytes_from_file(fileID, ifd_offset + 2, 12 * num_tags + 4);
            entries = containers.Map();
            
            for i = 1:num_tags
                entry = ifd_bytes(12 * (i - 1) + 1 : 12 * i);
                if verbose
                    fprintf("	Entry #%d:", i);
                    for e = 1:size(entry, 1)
                        fprintf(" %d,", entry(e));
                    end
                    fprintf("\n");
                end
                
                [tag_id, value] = ImageJ_formatted_TIFF.interpret_ifd_entry(fileID, entry, header.endian, verbose);
                if ~isempty(value)
                    entries(tag_id) = value;
                end
                
                try
                    eval(sprintf("header.%s = entries('%s');", tag_id, tag_id));
                    if tag_id == "ImageDescription"
                        header.NumberCharsInImageDescription = size(entries('ImageDescription'), 1);
                        header.OffsetOfImageDescription = ImageJ_formatted_TIFF.bytes_to_int(entry(9:12), header.endian);
                        header.ImageDescription = convertCharsToStrings(char(entries('ImageDescription')));
                    end
                    
                    if tag_id == "XResolution"
                        resolution_numerator = ImageJ_formatted_TIFF.bytes_to_int(header.XResolution(1:4), header.endian);
                        resolution_denominator = ImageJ_formatted_TIFF.bytes_to_int(header.XResolution(5:8), header.endian);                        
                        header.resolution =  round(single(resolution_denominator) / single(resolution_numerator), 5);
                    end
                catch
                    if verbose
                        warning("Simple_IFD class does not contain this Tiff tag: %d\n", tag_id);
                    end
                end
            end
        end
        
        function [tag_id, value] = interpret_ifd_entry(fileID, entry, endian, verbose)
            tag_id = ImageJ_formatted_TIFF.bytes_to_int(entry(1:2), endian);
            keySet = {254, 256, 257, 258, 259, 262, 270, 273, 277, 278, 279, 282, 283, 296, 339};
            valueSet = {'NewSubFileType', 'ImageWidth', 'ImageLength', 'BitsPerSample', 'Compression', 'PhotometricInterpretation', 'ImageDescription', 'StripOffsets', 'SamplesPerPixel', 'RowsPerStrip', 'StripByteCounts', 'XResolution', 'YResolution', 'ResolutionUnit', 'SampleFormat'};
            tag_id_lookup = containers.Map(keySet, valueSet);
            try
                tag_id = tag_id_lookup(tag_id);                
                
                data_type = ImageJ_formatted_TIFF.bytes_to_int(entry(3:4), endian);
                keySet = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
                valueSet = {{'BYTE', 1}, {'ASCII', 1}, {'SHORT', 2}, {'LONG', 4}, {'RATIONAL', 8}, {'SBYTE', 1}, {'UNDEFINED', 8}, {'SSHORT', 2}, {'SLONG', 4}, {'SRATIONAL', 8}, {'FLOAT', 4}, {'DOUBLE', 8}};
                data_type_lookup = containers.Map(keySet, valueSet);
                try
                    value = data_type_lookup(data_type);
                    data_type = string(value(1));
                    bytes_per_count = cell2mat(value(2));
                catch ME
                    warning("Unknown data type in TIF tag: %s", data_type);
                    throw(ME);
                end
                
                data_count = ImageJ_formatted_TIFF.bytes_to_int(entry(5:8), endian);
                value_size_bytes = data_count * bytes_per_count;
                if value_size_bytes <= 4
                    % The DataOffset directly encode the value
                    value = entry(9 : 8 + value_size_bytes);
                else
                    % The DataOffset encodes a pointer to the value
                    offset = ImageJ_formatted_TIFF.bytes_to_int(entry(9:12), endian);
                    value = ImageJ_formatted_TIFF.get_bytes_from_file(fileID, offset, value_size_bytes);
                end
                
                % We still haven't converted the value from bytes yet, but at least we
                % got the correct bytes that encode the value.
                switch data_type
                    case {'BYTE', 'SHORT', 'LONG'}
                        if data_count == 1
                            value = ImageJ_formatted_TIFF.bytes_to_int(value, endian);
                            if verbose
                                fprintf("	%s: %d\n", tag_id, value);
                            end
                        end
                    case 'ASCII'
                        if verbose
                            content = convertCharsToStrings(char(value));
                            fprintf("	%s:\n%s", tag_id, content);
                        end
                    case 'RATIONAL'
                        if verbose
                            fprintf("	%s:", tag_id);
                            for v = 1:size(value, 1)
                                fprintf(" %d,", value(v));
                            end
                            fprintf("\n");
                        end
                end
            catch
                if verbose
                    warning("Unknown tag ID in TIF: %d with value / offset: %d\n", tag_id, ImageJ_formatted_TIFF.bytes_to_int(entry(9:12), endian));
                end
                value = [];
            end
        end
        
        function bytes = get_bytes_from_file(fileID, offset, num_bytes)
            fseek(fileID, offset, 'bof');
            bytes = fread(fileID, num_bytes, 'uint8');
        end
        
        function value = bytes_to_int(bytes, endian)
            value = 0;
            if endian == 'L'
                for i = 1:size(bytes, 1)
                    value = value + bytes(i) * 256 ^ (i - 1);
                end
            elseif endian == 'B'
                for i = 1:size(bytes, 1)
                    value = value + bytes(i) * 256 ^ (size(bytes, 1) - i);
                end
            else
                ME = MException("'endian' must be either big or little");
                throw(ME)
            end
        end
        
    end
end

%% Test of write image stack
%     stack_in = uint16(linspace(0, 65535, 1920 * 1600 * 500));
%     stack_in = reshape(stack_in, [1920, 1600, 500]);
%     stack_in = permute(stack_in, [2 1 3]);
%     tic
%     fileID = fopen('test.tif', 'w+');
%     ImageJ_formatted_TIFF.write_IFD(fileID, 1920, 1600, class(stack_in), 1, 500, 1, 0.04, 39);
%     fwrite(fileID, matrix_in, class(stack_in));
%
%     fclose(fileID);
%     toc

%% Test of read image stack
%     fileID = fopen('test.tif', 'r');
%     header = ImageJ_formatted_TIFF.parse_tif(fileID);
%     tic
%     fseek(fileID, header.StripOffsets, 'bof');
%     matrix_out = fread(fileID, [1, 1920 * 1600], 'uint16');
%     fclose(fileID);
%     toc
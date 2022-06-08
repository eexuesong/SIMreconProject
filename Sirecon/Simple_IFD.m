classdef Simple_IFD
    %%
    %     An Image File Directory (IFD) is a collection of information similar to a header, and it is used to describe the bitmapped data to which it is attached.
    %     Like a header, it contains information on the height, width, and depth of the image, the number of color planes, and the type of data compression used on the bitmapped data.
    %
    %     One of the misconceptions about TIFF is that the information stored in the IFD tags is actually part of the TIFF header IFH.
    %     In fact, this information is often referred to as the "TIFF Header Information."
    %     While other formats do store the type of information found in the IFD in the header, the TIFF header does not contain this information.
    %     Therefore, it is possible to think of the IFDs in a TIFF file as extensions of the TIFF file header.
    %
    %     A TIFF file may contain any number of images, from zero on up:
    %       Each image is considered to be a separate subfile (i.e., a bitmap) and has an IFD describing the bitmapped data.
    %       Each TIFF subfile can be written as a separate TIFF file or can be stored with other subfiles in a single TIFF file.
    %       Each subfile bitmap and IFD may reside anywhere in the TIFF file after the headers, and there may be only one IFD per image.
    %     The last field of every IFD contains an offset value to the next IFD, if any.
    %     If the offset value of any IFD is 00h, then there are no more images left to read in the TIFF file.
    %
    %     The format of an Image File Directory (IFD) is shown in the following structure:
    %     typedef struct _TifIfd
    %     {
    %           WORD NumDirEntries; /* Number of Tags in IFD */
    %           TIFTAG TagList[]; /* Array of Tags */
    %           DWORD NextIFDOffset; /* Offset to next IFD */
    %     } TIFIFD;
    %
    %     NumDirEntries: is a 2-byte value indicating the number of tags found in the IFD.
    %                    Following this field is a series of tags; the number of tags corresponds to the value of the NumDirEntries field.
    %     TagList: Each tag structure is 12 bytes in size and, in the sample code above,
    %              is represented by an array of structures of the data type denition TIFTAG.
    %              The number of tags per IFD is limited to 65,535.
    %     NextIFDOffset: contains the offset position of the beginning of the next IFD.
    %                    If there are no more IFDs, then the value of this field is 00h.
    %
    %
    %     Tags: a tag can be thought of as a data field in a file header.
    %     Difference: A header field used to hold a byte of data need only be a byte in size.
    %                 A tag containing one byte of information, however, must always be twelve bytes in size.
    %     A TIFF tag has the following 12-byte structure:
    %     typedef struct _TifTag
    %     {
    %           WORD TagId; /* The tag identifier, 2 bytes */
    %           WORD DataType; /* The scalar type of the data items, 2 bytes */
    %           DWORD DataCount; /* The number of items in the tag data, 4 bytes */
    %           DWORD DataOffset; /* The byte offset to the data items, 4 bytes */
    %     } TIFTAG;
    %
    %     TagId: is a numeric value identifying the type of information the tag contains.
    %            More specically, the TagId indicates what the tag information represents.
    %            Typical information found in every TIFF file includes:
    %                   --- the height and width of the image, --- the depth of each pixel, and --- the type of data encoding used to compress the bitmap.
    %            Tags are normally identied by their TagId value and should always be written to an IFD in ascending order of the values found in the TagId field.
    %
    %     DataType: contains a value indicating the scalar data type of the information found in the tag.
    %               The following values are supported:
    %                       1 BYTE      8-bit unsigned integer;
    %                       2 ASCII     8-bit, NULL-terminated string;
    %                       3 SHORT     16-bit unsigned integer;
    %                       4 LONG      32-bit unsigned integer;
    %                       5 RATIONAL  Two 32-bit unsigned integers.
    %               The BYTE, SHORT, and LONG data types correspond to the BYTE, WORD, and DWORD data types used.
    %               The ASCII data type contains strings of 7-bit ASCII character data, which are always NULL-terminated and may be padded out to an even length if necessary.
    %               The RATIONAL data type is actually two LONG values and is used to store the two components of a fractional value.
    %                   The 1st value stores the numerator, and the 2nd value stores the denominator.
    %
    %     DataCount: indicates the number of items referenced by the tag and doesn't show the actual size of the data itself.
    %                Therefore, a DataCount of 08h does not necessarily indicate that eight bytes of data exist in the tag.
    %                This value indicates that eight items exist for the data type specied by this tag.
    %                For example, a DataCount value of 08h and a DataType of 03h indicate that the tag data is eight contiguous 16-bit unsigned integers, a total of 16 bytes in size.
    %                             A DataCount of 28h and a DataType of 02h indicate an ASCII character string 40 bytes in length,
    %                                       including the NULL-terminator character, but not any padding if present.
    %                             And a DataCount of 01h and a DataType of 05h indicate a single RATIONAL value a total of 8 bytes in size.
    %
    %     DataOffset: is a 4-byte field that contains the offset location of the actual tag data within the TIFF file.
    %                 1. If the tag data is 4 bytes or less in size, the data may be found in this field.
    %                 2. If the tag data is greater than 4 bytes in size, then this field contains an offset to the position of the data in the TIFF file.
    %                 Packing data within the DataOffset field is an optimization within the TIFF specication and is not required to be performed.
    %                 Most data is typically stored outside the tag, occurring before or after the IFD.
    %
    %
    %
    %     Note that several of these tags have default values that are used if the tag does not actually appear in a TIFF file.
    %     Bi-level (formerly Class B) and Gray-scale (formerly Class G) TIFF files must contain the 13 tags listed.
    %     These tags must appear in all revision 5.0 and 6.0 TIFF files regardless of the type of image data stored.
    %     Minimum Required Tags for TIFF Class B and Class G:
    %     Tag Type      TagId (2 bytes)       Tag Name                      DataType               DataCount N           default
    %       254              00feh            NewSubfileTyp            dword (LONG, 4 bytes)            1                   0
    %                   Currently defined values for the bitmap are:        0 - Image is reduced of another TIFF image in this file
    %                                                                       1 - Image is a single page of a multi-page
    %                                                                       2 - Image is a transparency mask for another image in this file
    %
    %       256              0100h            ImageWidth               word or dword                    1               No default
    %                   The image's width, in pixels (X:horizontal). The number of columns in the image.
    %
    %       257              0101h            ImageLength              word or dword                    1               No default
    %                   The image's length (height) in pixels (Y:vertical). The number of rows in the image.
    %
    %       258              0102h            BitsPerSample            word (SHORT, 2 bytes)     SamplesPerPixel            1
    %                   Number of bits per sample (bit depth). Note that this tag allows a different number of bits per sample for each sample corresponding to a pixel.
    %                   For example, RGB color data could use a different number of bits per sample for each of the three color planes.
    %
    %       259              0103h            Compression                    word                       1                   1
    %                   1 = No compression, but pack data into bytes as tightly as possible, with no unused bits except at the end of a row.
    %                       The byte ordering of data >8 bits must be consistent with that specified in the TIFF file header (bytes 0 and 1).
    %                       Rows are required to begin on byte boundaries.
    %                   2 = CCITT Group 3 1-Dimensional Modified Huffman run length encoding.
    %                   3 = Facsimile-compatible CCITT  Group 3.
    %                   4 = Facsimile-compatible CCITT  Group 4.
    %                   5 = LZW Compression, for grayscale, mapped color, and full color images.
    %                   32773 = PackBits compression, a simple byte oriented run length scheme for 1-bit images.
    %
    %       262              0106h            PhotometricInterpretation      word                       1               No default
    %                   0 = For bilevel and grayscale images: 0 is imaged as white. 2**BitsPerSample-1 is imaged as black.
    %                       If GrayResponseCurve exists, it overrides the PhotometricInterpretation value.
    %                   1 = For bilevel and grayscale images: 0 is imaged as black. 2**BitsPerSample-1 is imaged as white.
    %                       If GrayResponseCurve exists, it overrides the PhotometricInterpretation value.
    %                   2 = RGB. In the RGB model, a color is described as a combination of the three primary colors of light (red, green, and blue) in particular concentrations.
    %                       For each of the three samples,  0 represents minimum intensity, and 2**BitsPerSample - 1 represents maximum intensity.
    %                       For PlanarConfiguration = 1, the samples are stored in the indicated order: first Red, then Green, then Blue.
    %                       For PlanarConfiguration = 2, the StripOffsets for the sample planes are stored in the indicated order:
    %                                                    first the Red sample plane StripOffsets, then the Green plane StripOffsets, then the Blue plane StripOffsets.
    %                   3 = "Palette color." In this mode, a color is described with a single sample.
    %                       The sample is used as an index into ColorMap. The sample is used to index into each of the red, green and blue curve tables to retrieve an RGB triplet defining an actual color.
    %                       When this PhotometricInterpretation value is used, the color response curves must also be supplied. SamplesPerPixel must be 1.
    %                   4 = Transparency Mask. This means that the image is used to define an irregularly shaped region of another image in the same TIFF file.
    %                       SamplesPerPixel and BitsPerSample must be 1. PackBits compression is recommended. The 1-bits define the interior of the region;
    %                       the 0-bits define the exterior of the region. The Transparency Mask must have the same ImageLength and ImageWidth as the main image.
    %
    %       273              0111h            StripOffsets             word or dword    = StripsPerImage for PlanarConfiguration equal to 1;   No default
    %                                                                                   = SamplesPerPixel * StripsPerImage for PlanarConfiguration equal to 2.
    %                   For each strip, the byte offset of that strip. The offset is specified with respect to the beginning of the TIFF file.
    %                   Note that this implies that each strip has a location independent of the locations of other strips. This feature may be useful for editing applications.
    %                   This field is the only way for a reader to find the image data, and hence must exist.
    %
    %       277              0115h            SamplesPerPixel                word                       1                   1
    %                   The number of samples per pixel. SamplesPerPixel is 1 for bilevel, grayscale, and palette color images.
    %                                                    SamplesPerPixel is 3 for RGB images.
    %
    %       278              0116h            RowsPerStrip             word or dword                    1               2**32 - 1 (effectively infinity. That is, the entire image is one strip. Recomended is a strip size of 8K.)
    %                   The number of rows per strip.  The image data is organized into strips for fast access to individual rows when the data is compressed
    %                                                - though this field is valid even if the data is not compressed.
    %                   Noted be Xuesong 2019/04/16: the RowsPerStrip value species the maximum value, and not the required value, of the number of rows per strip.
    %                                                Many TIFF files, in fact, store a single strip of data and specify an arbitrarily large RowsPerStrip value.
    %
    %       279              0117h            StripByteCounts          word or dword     = StripsPerImage for PlanarConfiguration equal to 1;   No default
    %                                                                                    = SamplesPerPixel * StripsPerImage for PlanarConfiguration equal to 2.
    %                   For each strip, the number of bytes in that strip. The existence of this field greatly
    %                                   simplifies the chore of buffering compressed data, if the strip size is reasonable.
    %
    %       282              011Ah            XResolution                  RATIONAL                     1               No default
    %                   The number of pixels per ResolutionUnit in the X direction, i.e., in the  ImageWidth direction
    %       283              011Bh            YResolution                  RATIONAL                     1               No default
    %                   The number of pixels per ResolutionUnit in the Y direction, i.e., in the ImageLength direction.
    %       296              0128h            ResolutionUnit                 word                       1                    2
    %                   1 = No absolute unit of measurement. Used for images that may have a non-square aspect ratio, but no meaningful absolute dimensions.
    %                       The drawback of ResolutionUnit = 1 is that different applications will import the image at different sizes.
    %                       Even if the decision is quite arbitrary, it might be better to use dots per inch or dots per centimeter,
    %                       and pick XResolution and YResolution such that the aspect ratio is correct and the maximum dimension of the image is about four inches.
    %                   2 = Inch.
    %                   3 = Centimeter.
    %
    %     Palette-color (formerly Class P) TIFF files add a 14th required tag that describes the type of palette information found within the TIFF image file:
    %       320              0140h            ColorMap
    %
    %     RGB (formerly Class R) TIFF files contain the same tags as bi-level TIFF (Class B) files and add a 14th required tag, which describes the format of the bitmapped data in the image:
    %       284              011ch            PlanarConfiguration            word                       1                    1
    %                   1 = The sample values for each pixel are stored contiguously, so that there is a single image plane.
    %                       See PhotometricInterpretation to determine the order of the samples within the pixel data.
    %                       So, for RGB data, the data is stored RGBRGBRGB...and so on.
    %                   2 = The samples are stored in separate "sample planes."
    %                       The values in StripOffsets and StripByteCounts are then arranged as a 2-dimensional array, with SamplesPerPixel rows and StripsPerImage columns.
    %                       (All of the columns for row 0 are stored first, followed by the columns of row 1, and so on.)
    %                       PhotometricInterpretation describes the type of data that is stored in each sample plane.
    %                       For example,  RGB data is stored with the Red samples in one sample plane, the Green in another, and the Blue in another.
    %                   Noted by Xuesong 2019/04/16: If SamplesPerPixel is 1 (bilevel, grayscale),
    %                                                PlanarConfiguration is irrelevant, and should not be included.
    %
    %     YCbCr TIFF files add 4 additional tags to the baseline:
    %       529              0211h            YCbCrCoecients
    %       530              0212h            YCbCrSubSampling
    %       531              0213h            YCbCrPositioning
    %       532              0214h            ReferenceBlackWhite
    %
    %
    %     In Andy's code, he also used the following tags:
    %       270              010eh            ImageDescription              ASCII
    %                   For example, a user may wish to attach a comment such as "1988 company picnic" to an image.
    %       339              0153h            SampleFormat             word (SHORT, 2 bytes)     SamplesPerPixel   1 (unsigned integer data)
    %                   Specifies how to interpret each data sample in a pixel. The specification defines these values:
    %                   1 = unsigned integer data;                  SAMPLEFORMAT_UINT = 1;
    %                   2 = two's complement signed integer data;   SAMPLEFORMAT_INT = 2;
    %                   3 = IEEE floating point data;               SAMPLEFORMAT_IEEEFP = 3;
    %                   4 = undefined data format.                  SAMPLEFORMAT_VOID = 4;
    %                   Noted by Xuesong 2019/04/16: SampleFormat field does not specify the size of data samples;
    %                                                this is still done by the BitsPerSample field.
    %
    %     And he did not define the following required tags:
    %       259              0103h            Compression           Default = 1, No compression
    %       282              011Ah            XResolution           No default;
    %       283              011Bh            YResolution           No default; That is probably the reason why iSIM tiff data cannot be directed read by Windows. By Xuesong 2019/04/16
    %       296              0128h            ResolutionUnit        Default = 2, Inch
    
    
    %% A very simple TIF IFD with 15 tags (2 + 15*12 + 4 = 186 bytes)
    properties(Dependent)
        bytes
    end
    
    properties
        NumDirEntries = 15
        
        NewSubFileType = 0
        ImageWidth = 0
        ImageLength = 0
        BitsPerSample = 0
        Compression = 1
        PhotometricInterpretation = 1
        NumberCharsInImageDescription = 0
        OffsetOfImageDescription = 0
        StripOffsets = 0
        SamplesPerPixel = 1
        RowsPerStrip = 0
        StripByteCounts = 0
        XResolution = 0
        YResolution = 0
        ResolutionUnit = 3
        SampleFormat = 1
        NextIFD = 0
        
        endian = "L"
        ImageDescription
        resolution
    end
    
    methods
        function obj = Simple_IFD(endian)
            if nargin == 0
                obj.endian = 'L';
            else
                obj.endian = endian;
            end
        end
        
        function bytes=get.bytes(obj)
            % A very simple TIF IFD with 15 tags (2 + 15*12 + 4 = 186 bytes)
            [~, ~, system_endian] = computer;
            if system_endian == obj.endian
                swapbytes_flag = 0;
            else
                swapbytes_flag = 1;
            end
            
            bytes = zeros([1 2], 'uint8');
            if swapbytes_flag == 0
                bytes(1:2) = typecast(uint16(obj.NumDirEntries), 'uint8');
            else
                bytes(1:2) = typecast(swapbytes(uint16(obj.NumDirEntries)), 'uint8');
            end
            
            %% NewSubFileType: TagId = 254 (FE, 00); DataType = 4 (LONG); DataCount = 1; DataOffset = 0
            % Little endian:    254,   0,   4,   0,   1,   0,   0,   0,   0,   0,   0,   0
            % Big endian:       0,   254,   0,   4,   0,   0,   0,   1,   0,   0,   0,   0
            TagList = Simple_IFD.tag_array(254, 4, 1, obj.NewSubFileType, swapbytes_flag);
            bytes = [bytes, TagList];
            
            %% ImageWidth: TagId = 256 (00, 01); DataType = 4 (LONG); DataCount = 1; DataOffset = 0
            % Little endian:    0,   1,   4,   0,   1,   0,   0,   0,   0,   0,   0,   0
            % Big endian:       1,   0,   0,   4,   0,   0,   0,   1,   0,   0,   0,   0
            TagList = Simple_IFD.tag_array(256, 4, 1, obj.ImageWidth, swapbytes_flag);
            bytes = [bytes, TagList];
            
            %% ImageLength: TagId = 257 (01, 01); DataType = 4 (LONG); DataCount = 1; DataOffset = 0
            % Little endian:    1,   1,   4,   0,   1,   0,   0,   0,   0,   0,   0,   0
            % Big endian:       1,   1,   0,   4,   0,   0,   0,   1,   0,   0,   0,   0
            TagList = Simple_IFD.tag_array(257, 4, 1, obj.ImageLength, swapbytes_flag);
            bytes = [bytes, TagList];
            
            %% BitsPerSample: TagId = 258 (02; 01); DataType = 3 (SHORT); DataCount = 1; DataOffset = 0
            % Little endian:    2,   1,   3,   0,   1,   0,   0,   0,   0,   0,   0,   0
            % Big endian:       1,   2,   0,   3,   0,   0,   0,   1,   0,   0,   0,   0
            TagList = Simple_IFD.tag_array(258, 3, 1, obj.BitsPerSample, swapbytes_flag);
            bytes = [bytes, TagList];
            
            %% Compression: TagId = 259 (03; 01); DataType = 3 (SHORT); DataCount = 1; DataOffset = 1 (No compression)
            % Little endian:    3,   1,   3,   0,   1,   0,   0,   0,   1,   0,   0,   0
            % Big endian:       1,   3,   0,   3,   0,   0,   0,   1,   0,   1,   0,   0
            TagList = Simple_IFD.tag_array(259, 3, 1, obj.Compression, swapbytes_flag);
            bytes = [bytes, TagList];
            
            %% PhotometricInterpretation: TagId = 262 (06; 01); DataType = 3 (SHORT); DataCount = 1; DataOffset = 1 (0 is imaged as black. 2**BitsPerSample-1 is imaged as white)
            % Little endian:    6,   1,   3,   0,   1,   0,   0,   0,   1,   0,   0,   0
            % Big endian:       1,   6,   0,   3,   0,   0,   0,   1,   0,   1,   0,   0
            TagList = Simple_IFD.tag_array(262, 3, 1, obj.PhotometricInterpretation, swapbytes_flag);
            bytes = [bytes, TagList];
            
            %% ImageDescription: TagId = 270 (0E; 01); DataType = 2 (ASCII); DataCount = 0; DataOffset = 0
            % Little endian:    14,   1,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0
            % Big endian:       1,   14,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0
            TagList = Simple_IFD.tag_array(270, 2, obj.NumberCharsInImageDescription, obj.OffsetOfImageDescription, swapbytes_flag);
            bytes = [bytes, TagList];
            
            %% StripOffsets: TagId = 273 (11; 01); DataType = 4 (LONG); DataCount = 1; DataOffset = 0
            % Little endian:    17,   1,   4,   0,   1,   0,   0,   0,   0,   0,   0,   0
            % Big endian:       1,   17,   0,   4,   0,   0,   0,   1,   0,   0,   0,   0
            TagList = Simple_IFD.tag_array(273, 4, 1, obj.StripOffsets, swapbytes_flag);
            bytes = [bytes, TagList];
            
            %% SamplesPerPixel: TagId = 277 (15; 01); DataType = 3 (SHORT); DataCount = 1; DataOffset = 1 (bilevel; grayscale; and palette color images)
            % Little endian:    21,   1,   3,   0,   1,   0,   0,   0,   1,   0,   0,   0
            % Big endian:       1,   21,   0,   3,   0,   0,   0,   1,   0,   1,   0,   0
            TagList = Simple_IFD.tag_array(277, 3, 1, obj.SamplesPerPixel, swapbytes_flag);
            bytes = [bytes, TagList];
            
            %% RowsPerStrip: TagId = 278 (16; 01); DataType = 3 (SHORT); DataCount = 1; DataOffset = 0
            % Little endian:    22,   1,   3,   0,   1,   0,   0,   0,   0,   0,   0,   0
            % Big endian:       1,   22,   0,   3,   0,   0,   0,   1,   0,   0,   0,   0
            TagList = Simple_IFD.tag_array(278, 3, 1, obj.RowsPerStrip, swapbytes_flag);
            bytes = [bytes, TagList];
            
            %% StripByteCounts: TagId = 279 (17; 01); DataType = 4 (LONG); DataCount = 1; DataOffset = 0
            % Little endian:    23,   1,   4,   0,   1,   0,   0,   0,   0,   0,   0,   0
            % Big endian:       1,   23,   0,   4,   0,   0,   0,   1,   0,   0,   0,   0
            TagList = Simple_IFD.tag_array(279, 4, 1, obj.StripByteCounts, swapbytes_flag);
            bytes = [bytes, TagList];
            
            %% XResolution: TagId = 282 (1A; 01); DataType = 5 (RATIONAL); DataCount = 1; DataOffset = 0
            % Little endian:    23,   1,   4,   0,   1,   0,   0,   0,   0,   0,   0,   0
            % Big endian:       1,   23,   0,   4,   0,   0,   0,   1,   0,   0,   0,   0
            TagList = Simple_IFD.tag_array(282, 5, 1, obj.XResolution, swapbytes_flag);
            bytes = [bytes, TagList];
            
            %% YResolution: TagId = 283 (1B; 01); DataType = 5 (RATIONAL); DataCount = 1; DataOffset = 0
            % Little endian:    27,   1,   5,   0,   1,   0,   0,   0,   0,   0,   0,   0
            % Big endian:       1,   27,   0,   5,   0,   0,   0,   1,   0,   0,   0,   0
            TagList = Simple_IFD.tag_array(283, 5, 1, obj.YResolution, swapbytes_flag);
            bytes = [bytes, TagList];
            
            %% ResolutionUnit: TagId = 296 (28; 01); DataType = 3 (SHORT); DataCount = 1; DataOffset = 3 (Centimeter)
            % Little endian:    40,   1,   3,   0,   1,   0,   0,   0,   3,   0,   0,   0
            % Big endian:       1,   40,   0,   3,   0,   0,   0,   1,   0,   3,   0,   0
            TagList = Simple_IFD.tag_array(296, 3, 1, obj.ResolutionUnit, swapbytes_flag);
            bytes = [bytes, TagList];
            
            %% SampleFormat: TagId = 339 (53; 01); DataType = 3 (SHORT); DataCount = 1; DataOffset = 0
            % Little endian:    83,   1,   3,   0,   1,   0,   0,   0,   0,   0,   0,   0
            % Big endian:       1,   83,   0,   3,   0,   0,   0,   1,   0,   0,   0,   0
            TagList = Simple_IFD.tag_array(339, 3, 1, obj.SampleFormat, swapbytes_flag);
            bytes = [bytes, TagList];
            
            NextIFDOffset = zeros([1 4], 'uint8');
            if swapbytes_flag == 0
                NextIFDOffset(1:4) = typecast(uint32(obj.NextIFD), 'uint8');
            else
                NextIFDOffset(1:4) = typecast(swapbytes(uint32(obj.NextIFD)), 'uint8');
            end
            
            bytes = [bytes, NextIFDOffset];
            
            %             bytes = transpose(uint8([15;    0;
            %             %NewSubFileType: TagId = 254 (FE, 00); DataType = 4 (LONG); DataCount = 1; DataOffset = 0
            %             254;   0;   4;   0;   1;   0;   0;   0;   0;   0;   0;   0;
            %             % ImageWidth: TagId = 256 (00; 01); DataType = 4 (LONG); DataCount = 1; DataOffset = 0
            %             0;   1;   4;   0;   1;   0;   0;   0;   0;   0;   0;   0;
            %             % ImageLength: TagId = 257 (01; 01); DataType = 4 (LONG); DataCount = 1; DataOffset = 0
            %             1;   1;   4;   0;   1;   0;   0;   0;   0;   0;   0;   0;
            %             % BitsPerSample: TagId = 258 (02; 01); DataType = 3 (SHORT); DataCount = 1; DataOffset = 0
            %             2;   1;   3;   0;   1;   0;   0;   0;   0;   0;   0;   0;
            %             % Compression: TagId = 259 (03; 01); DataType = 3 (SHORT); DataCount = 1; DataOffset = 1 (No compression)
            %             3;   1;   3;   0;   1;   0;   0;   0;   1;   0;   0;   0;
            %             % PhotometricInterpretation: TagId = 262 (06; 01); DataType = 3 (SHORT); DataCount = 1; DataOffset = 1 (0 is imaged as black. 2**BitsPerSample-1 is imaged as white)
            %             6;   1;   3;   0;   1;   0;   0;   0;   1;   0;   0;   0;
            %             % ImageDescription: TagId = 270 (0E; 01); DataType = 2 (ASCII); DataCount = 0; DataOffset = 0
            %             14;   1;   2;   0;   0;   0;   0;   0;   0;   0;   0;   0;
            %             % StripOffsets: TagId = 273 (11; 01); DataType = 4 (LONG); DataCount = 1; DataOffset = 0
            %             17;   1;   4;   0;   1;   0;   0;   0;   0;   0;   0;   0;
            %             % SamplesPerPixel: TagId = 277 (15; 01); DataType = 3 (SHORT); DataCount = 1; DataOffset = 1 (bilevel; grayscale; and palette color images)
            %             21;   1;   3;   0;   1;   0;   0;   0;   1;   0;   0;   0;
            %             % RowsPerStrip: TagId = 278 (16; 01); DataType = 3 (SHORT); DataCount = 1; DataOffset = 0
            %             22;   1;   3;   0;   1;   0;   0;   0;   0;   0;   0;   0;
            %             % StripByteCounts: TagId = 279 (17; 01); DataType = 4 (LONG); DataCount = 1; DataOffset = 0
            %             23;   1;   4;   0;   1;   0;   0;   0;   0;   0;   0;   0;
            %             % XResolution: TagId = 282 (1A; 01); DataType = 5 (RATIONAL); DataCount = 1; DataOffset = 0
            %             26;   1;   5;   0;   1;   0;   0;   0;   0;   0;   0;   0;
            %             % YResolution: TagId = 283 (1B; 01); DataType = 5 (RATIONAL); DataCount = 1; DataOffset = 0
            %             27;   1;   5;   0;   1;   0;   0;   0;   0;   0;   0;   0;
            %             % ResolutionUnit: TagId = 296 (28; 01); DataType = 3 (SHORT); DataCount = 1; DataOffset = 3 (Centimeter)
            %             40;   1;   3;   0;   1;   0;   0;   0;   3;   0;   0;   0;
            %             % SampleFormat: TagId = 339 (53; 01); DataType = 3 (SHORT); DataCount = 1; DataOffset = 0
            %             83;   1;   3;   0;   1;   0;   0;   0;   0;   0;   0;   0;
            %             % Next IFD = 0
            %             0;   0;   0;   0;]));
            %             bytes(11:14) = typecast(uint32(obj.NewSubFileType), 'uint8');
            %             bytes(23:26) = typecast(uint32(obj.ImageWidth), 'uint8');
            %             bytes(35:38) = typecast(uint32(obj.ImageLength), 'uint8');
            %             bytes(47:50) = typecast(uint32(obj.BitsPerSample), 'uint8');
            %             bytes(59:62) = typecast(uint32(obj.Compression), 'uint8');
            %             bytes(71:74) = typecast(uint32(obj.PhotometricInterpretation), 'uint8');
            %             bytes(79:82) = typecast(uint32(obj.NumberCharsInImageDescription), 'uint8');
            %             bytes(83:86) = typecast(uint32(obj.OffsetOfImageDescription), 'uint8');
            %             bytes(95:98) = typecast(uint32(obj.StripOffsets), 'uint8');
            %             bytes(107:110) = typecast(uint32(obj.SamplesPerPixel), 'uint8');
            %             bytes(119:122) = typecast(uint32(obj.RowsPerStrip), 'uint8');
            %             bytes(131:134) = typecast(uint32(obj.StripByteCounts), 'uint8');
            %             bytes(143:146) = typecast(uint32(obj.XResolution), 'uint8');
            %             bytes(155:158) = typecast(uint32(obj.YResolution), 'uint8');
            %             bytes(167:170) = typecast(uint32(obj.ResolutionUnit), 'uint8');
            %             bytes(179:182) = typecast(uint32(obj.SampleFormat), 'uint8');
            %             bytes(183:186) = typecast(uint32(obj.NextIFD), 'uint8');
        end
        
        function obj = set_dtype(obj, dtype)
            keySet = {'uint8', 'uint16', 'uint32', 'uint64', 'int8', 'int16', 'int32', 'int64', 'single', 'double'};
            valueSet = {[1, 8], [1, 16],  [1, 32], [1, 64], [2, 8], [2, 16], [2, 32], [2, 64], [3, 32], [3, 64]};
            allowed_dtypes = containers.Map(keySet, valueSet);
            try
                value = allowed_dtypes(dtype);
                obj.SampleFormat = value(1);
                obj.BitsPerSample = value(2);
            catch ME
                warning_msg = sprintf("Array datatype '%s' not allowed. Allowed types:\n", dtype);
                for i = 1:length(allowed_dtypes)
                    warning_msg = sprintf("%s%s\n", warning_msg, keySet{i});
                end
                warning(warning_msg);
                throw(ME);
            end
        end
    end
    
    methods(Static)
        function TagList = tag_array(TagId, DataType, DataCount, DataOffset, Swapbytes_Flag)
            keySet = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
            valueSet = {{'BYTE', 1}, {'ASCII', 1}, {'SHORT', 2}, {'LONG', 4}, {'RATIONAL', 8}, {'SBYTE', 1}, {'UNDEFINED', 8}, {'SSHORT', 2}, {'SLONG', 4}, {'SRATIONAL', 8}, {'FLOAT', 4}, {'DOUBLE', 8}};
            data_type_lookup = containers.Map(keySet, valueSet);
            try
                value = data_type_lookup(DataType);
                data_type = string(value(1));
                bytes_per_count = cell2mat(value(2));
            catch ME
                warning("Unknown data type in TIF tag: %d", DataType);
                throw(ME);
            end
            value_size_bytes = DataCount * bytes_per_count;
            
            TagList = zeros([1 12], 'uint8');
            if Swapbytes_Flag == 0
                TagList(1:2) = typecast(uint16(TagId), 'uint8');
                TagList(3:4) = typecast(uint16(DataType), 'uint8');
                TagList(5:8) = typecast(uint32(DataCount), 'uint8');
                switch value_size_bytes
                    case 1
                        TagList(9) = uint8(DataOffset);
                    case 2
                        TagList(9 : 10) = typecast(uint16(DataOffset), 'uint8');
                    otherwise
                        TagList(9 : 12) = typecast(uint32(DataOffset), 'uint8');
                end
            else
                TagList(1:2) = typecast(swapbytes(uint16(TagId)), 'uint8');
                TagList(3:4) = typecast(swapbytes(uint16(DataType)), 'uint8');
                TagList(5:8) = typecast(swapbytes(uint32(DataCount)), 'uint8');                
                switch value_size_bytes
                    case 1
                        TagList(9) = uint8(DataOffset);
                    case 2
                        TagList(9 : 10) = typecast(swapbytes(uint16(DataOffset)), 'uint8');
                    otherwise
                        TagList(9 : 12) = typecast(swapbytes(uint32(DataOffset)), 'uint8');
                end
            end
        end
    end
end


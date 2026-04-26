function w = decompressShape(shape, forceDecompression)
%decompressShape Decompress a gradient or pulse shape.
%
%   PURPOSE
%     Invert the run-length-on-derivative encoding produced by
%     mr.compressShape, returning the original sampled waveform as a
%     column vector. The input is the struct stored in the [SHAPES]
%     section of a Pulseq .seq file (or its in-memory equivalent inside
%     a shape library entry). Used internally by mr.Sequence when reading
%     a sequence file or reconstructing waveforms via getBlock, and
%     callable directly when inspecting a stored shape.
%
%   SIGNATURES
%     w = mr.decompressShape(shape)
%     w = mr.decompressShape(shape, forceDecompression)
%
%     By default the function detects an uncompressed payload by checking
%     length(shape.data) == shape.num_samples and returns shape.data as a
%     column vector untouched. Set forceDecompression to true to skip
%     that shortcut and always run the run-length decoder; this is
%     required when reading legacy .seq files in which a genuinely
%     compressed payload happens to satisfy the length-equality check.
%
%   INPUTS
%     shape              [required]    struct produced by mr.compressShape (or
%                                      read from a [SHAPES] block) with fields:
%                                        .num_samples : scalar, number of samples
%                                                       in the uncompressed waveform
%                                        .data        : row vector, either the
%                                                       compressed encoding or the
%                                                       original samples (see NOTES)
%     forceDecompression [optional]    logical scalar. If true, always run the
%                                      run-length decoder even when length(data)
%                                      equals num_samples. Default: false.
%
%   OUTPUT
%     w  column vector, double, length shape.num_samples. Same numeric units
%        as the waveform originally passed to mr.compressShape (gradient
%        amplitude in Hz/m, RF magnitude in Hz, RF phase in radians, or a
%        time vector in raster units, depending on what was compressed).
%
%   ERRORS
%     - 'MATLAB:nonExistentField': shape is missing the .data or
%       .num_samples field.
%     - assertion failure (no message): the compressed payload is
%       malformed and the trailing unpacked tail length does not match
%       num_samples - countUnpack. Indicates a corrupt .seq file or a
%       payload that should have been read with forceDecompression=true.
%
%   NOTES
%     - The output is always a column vector regardless of the input
%       orientation, because the function ends with w = w(:).
%     - The shortcut for uncompressed payloads is purely a length check.
%       In Pulseq v1.4.x and later the encoder guarantees that
%       length(data) == num_samples implies the payload is the raw
%       waveform, so the default forceDecompression=false is correct for
%       any sequence written by current Pulseq tools. Older files may
%       require forceDecompression=true to be decoded faithfully.
%     - Round-trip accuracy is bounded by the 1e-7 quantization step
%       used in mr.compressShape, not by floating-point epsilon. Do not
%       rely on bitwise equality between the input to compressShape and
%       the output of decompressShape.
%     - The function does not consult any system or raster information.
%       It returns dimensionless samples; converting to physical time
%       (e.g., multiplying an RF time-shape result by rfRasterTime) is
%       the caller's responsibility, as done in mr.Sequence/getBlock.
%
%   EXAMPLE
%     % Compress a trapezoidal gradient waveform, then decompress it and
%     % verify the round-trip.
%     sys = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
%                   'MaxSlew', 170, 'SlewUnit', 'T/m/s');
%     g = mr.makeTrapezoid('x', 'Area', 1000, 'system', sys);
%     raster = sys.gradRasterTime;
%     nRise = round(g.riseTime/raster);
%     nFlat = round(g.flatTime/raster);
%     nFall = round(g.fallTime/raster);
%     w = g.amplitude * [linspace(0,1,nRise+1), ...
%                        ones(1, nFlat-1), ...
%                        linspace(1,0,nFall+1)];
%     s = mr.compressShape(w(:));
%     w_recon = mr.decompressShape(s);
%     fprintf('samples: %d, max round-trip error: %.2e\n', ...
%             length(w_recon), max(abs(w_recon(:) - w(:))));
%
%   SEE ALSO
%     mr.compressShape, mr.pts2waveform


dataPack = shape.data;
dataPackLen = length(dataPack);
numSamples=shape.num_samples;

if nargin<2
    forceDecompression=false;
end

if ~forceDecompression && numSamples==dataPackLen
    % uncompressed shape
    w=dataPack';
    return;
end

w= zeros(1, numSamples) ;                 % pre-allocate the result matrix
                                          % dimensons: (1,length of the data set)
                                       
% decompression starts here

dataPackDiff = dataPack(2:end) - dataPack(1:end-1);

% when dataPackDiff == 0 the subsequent samples are equal ==> marker for
% repeats (run-length encoding)
dataPackMarkers=find(dataPackDiff==0.0);

countPack= 1;                                               % counter 1: points to the current compressed sample
countUnpack= 1;                                             % counter 2: points to the current uncompressed sample

for i=1:length(dataPackMarkers)
    nextPack=dataPackMarkers(i); % careful, this index may have "false positives" , e.g. if the value 3 repeats 3 times, then we will have 3 3 3
    currUnpackSamples=nextPack-countPack;
    if currUnpackSamples < 0 % this rejects false positives
        continue;        
    elseif currUnpackSamples > 0 % do we have an unpacked block to copy?
        w(countUnpack:(countUnpack+currUnpackSamples-1)) = dataPack(countPack:(nextPack-1));
        countPack = countPack + currUnpackSamples;
        countUnpack = countUnpack + currUnpackSamples;
    end
    % now comes the packed/repeated section
    rep=dataPack(countPack+2)+2;
    w(countUnpack:(countUnpack+rep-1))=dataPack(countPack);
    countPack= countPack + 3;
    countUnpack= countUnpack + rep;
end

% samples left?
if (countPack<=dataPackLen)
    assert(dataPackLen-countPack==numSamples-countUnpack);
    % copy the rest of the shape, it is unpacked
    w(countUnpack:end)= dataPack(countPack:end);
end

w = cumsum(w);
w=w(:);

% % test the new function against the old slow version
% w1=mr.decompressShape0(shape);
% assert(all(w==w1));

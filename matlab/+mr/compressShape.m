function s=compressShape(w, forceCompression)
%compressShape Compress a gradient or pulse shape.
%
%   PURPOSE
%     Encode a sampled waveform (gradient amplitude, RF magnitude, RF
%     phase, or a time vector in raster units) into the run-length
%     derivative format used by Pulseq .seq files. The scheme stores the
%     first sample followed by quantized derivatives, and collapses runs
%     of three or more equal derivatives into a (value, value, repeats)
%     triplet. This makes constant and piecewise-linear waveforms (the
%     common gradient and time-shape cases) extremely compact while
%     leaving arbitrary waveforms only marginally larger than the input.
%     The returned struct is consumed by mr.Sequence/write when emitting
%     the [SHAPES] section, and inverted by mr.decompressShape on read.
%
%   SIGNATURES
%     s = mr.compressShape(w)
%     s = mr.compressShape(w, forceCompression)
%
%     If length(w) <= 4 and forceCompression is false (the default), the
%     waveform is stored uncompressed (s.data = w). Set forceCompression
%     to true to apply the run-length encoding regardless of length.
%     For longer waveforms the compressed encoding is used only if it is
%     shorter than the original; otherwise the original samples are kept.
%
%   INPUTS
%     w                 [required]    numeric vector, the waveform samples to
%                                     compress (any unit; treated as a plain
%                                     numeric sequence). Must be finite. May be
%                                     a row or column vector; a column vector is
%                                     internally flattened to a row.
%     forceCompression  [optional]    logical scalar. If true, always emit the
%                                     run-length encoding even for short or
%                                     incompressible inputs. Default: false.
%
%   OUTPUT
%     s  struct with fields (in order returned by fieldnames):
%       .num_samples  scalar double, the number of samples in the original
%                     uncompressed waveform (i.e., numel(w)).
%       .data         row vector, double. Either the compressed encoding
%                     (first sample, then quantized derivatives with run
%                     lengths) or, when compression is skipped or would
%                     not shrink the data, the original samples as a row
%                     vector. mr.decompressShape distinguishes the two
%                     cases by comparing length(s.data) with s.num_samples.
%
%   ERRORS
%     - 'compressShape() received infinite samples': any element of w is
%       not finite (Inf or NaN).
%
%   NOTES
%     - Quantization step is 1e-7 (about 7 decimal digits, matching the
%       precision of single-precision floats and of the Pulseq .seq text
%       format). Round-trip error through mr.decompressShape is bounded
%       by this quantum, not by floating-point epsilon. Do not rely on
%       bitwise round-trip equality.
%     - Values of magnitude smaller than 1e-10 in the compressed stream
%       are flushed to exactly 0 to avoid spurious tiny derivatives.
%     - The "is compression worthwhile" test compares lengths of the
%       compressed and original streams. For mostly-random inputs the
%       function silently falls back to storing the original samples;
%       inspect length(s.data) versus s.num_samples to tell which branch
%       was taken (equal => uncompressed; otherwise compressed).
%     - Pass forceCompression = true when the consumer requires the
%       compressed format unconditionally (e.g., to keep all shape
%       entries on a single decoding path).
%
%   EXAMPLE
%     % Sample a trapezoidal gradient on the gradient raster, compress it
%     % into Pulseq shape format, and verify the round-trip.
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
%     fprintf('original samples: %d, stored samples: %d\n', ...
%             s.num_samples, length(s.data));
%     w_recon = mr.decompressShape(s);
%     fprintf('max round-trip error: %.2e\n', max(abs(w_recon(:) - w(:))));
%
%   SEE ALSO
%     mr.decompressShape, mr.pts2waveform

if nargin<2
    forceCompression=false;
end

if any(~isfinite(w))
    error('compressShape() received infinite samples');
end

if ~forceCompression && length(w) <= 4 % avoid compressing very short shapes
    s.num_samples=length(w);
    s.data = w(:)';
    return;
end


% %MZ: old code with implicit quantization
% data = [w(1); diff(w(:))];
% maskChanges = [true; abs(diff(data))>1e-8];   % TRUE if values change
% vals = data(maskChanges);                     % Elements without repetitions

% MZ: explicit quantization with error correction
quant_fac=1e-7; % single precision floating point has ~7.25 decimal places 
ws=w./quant_fac;
datq=round([ws(1); diff(ws(:))]);
qerr=ws(:)-cumsum(datq);
qcor=[0; diff(round(qerr))];
datd=datq+qcor;
maskChanges=[true; diff(datd)~=0];
vals=datd(maskChanges).*quant_fac;            % Elements without repetitions

k = find([maskChanges', true]);               % Indices of changes
n = diff(k)';                                 % Number of repetitions

% Encode in Pulseq format
nExtra=n-2;
vals2=vals; 
vals2(nExtra<0)=nan;
nExtra(nExtra<0)=nan;
v=[vals vals2 nExtra]';
v=v(isfinite(v));
v(abs(v)<1e-10)=0;
s.num_samples = length(w);
% decide whether compression makes sense, otherwise store the original
if forceCompression || s.num_samples > length(v)
    s.data = v';
else
    s.data = w(:)';
end



function grad=makeArbitraryGrad(channel,varargin)
%makeArbitraryGrad Create a gradient event with an arbitrary waveform.
%
%   PURPOSE
%     Build a gradient event from a user-supplied waveform sampled on the
%     gradient raster. The waveform is checked against gradient amplitude
%     and slew-rate hardware limits, and the boundary samples ('first' and
%     'last' values, just outside the shape) are recorded so the sequence
%     can connect this gradient seamlessly to the surrounding blocks. The
%     returned struct is consumed by mr.Sequence/addBlock to add the
%     gradient to a sequence.
%
%   SIGNATURES
%     g = mr.makeArbitraryGrad(channel, waveform)
%     g = mr.makeArbitraryGrad(channel, waveform, system)
%     g = mr.makeArbitraryGrad(channel, waveform, 'system', system, ...)
%     g = mr.makeArbitraryGrad(channel, waveform, system, 'first', f0, 'last', f1, ...)
%
%     The waveform is interpreted as samples on the gradient raster
%     (system.gradRasterTime). When 'oversampling' is true the waveform
%     is on the 1/2-raster grid and must have an odd number of samples.
%     Always supply 'first' and 'last' explicitly: omitting them triggers
%     a deprecation warning and the boundary values are extrapolated from
%     the first/last two samples of the waveform.
%
%   INPUTS
%     channel       [required]    char, gradient axis: 'x', 'y', or 'z'.
%     waveform      [required]    numeric vector, gradient amplitude in Hz/m,
%                                 one sample per gradient raster (or per
%                                 1/2-raster if oversampling=true). Reshaped
%                                 internally to a column vector.
%     system        [optional]    struct from mr.opts. May be passed positionally
%                                 OR as 'system', value. If omitted, mr.opts() is
%                                 called for defaults.
%     'oversampling' [name/value] logical, default false. If true the waveform
%                                 is sampled at 1/2 the gradient raster time
%                                 (length must be odd).
%     'maxGrad'     [name/value]  numeric, Hz/m, default 0 (use system.maxGrad).
%                                 Override gradient-amplitude limit for this call.
%     'maxSlew'     [name/value]  numeric, Hz/m/s, default 0 (use system.maxSlew).
%                                 Override slew-rate limit for this call.
%     'delay'       [name/value]  numeric, seconds, default 0. Delay from the
%                                 start of the block to the first waveform sample.
%     'first'       [name/value]  numeric, Hz/m, default NaN. Gradient amplitude
%                                 at the leading edge of the shape (one half
%                                 raster before the first sample). NaN triggers
%                                 the extrapolation warning.
%     'last'        [name/value]  numeric, Hz/m, default NaN. Gradient amplitude
%                                 at the trailing edge of the shape (one half
%                                 raster after the last sample). NaN triggers
%                                 the extrapolation warning.
%
%   OUTPUT
%     grad struct with fields (in order returned by fieldnames):
%       .type       char, always 'grad'
%       .channel    char, 'x' | 'y' | 'z'
%       .waveform   column vector, Hz/m, the gradient samples
%       .delay      double, seconds, leading delay before the first sample
%       .area       double, 1/m, integral of the waveform over its duration
%                   (sum(waveform).*gradRasterTime; halved-by-decimation when
%                   oversampling is active so the value reflects the physical
%                   gradient moment)
%       .tt         column vector, seconds, sample time stamps relative to
%                   the start of the shape. Without oversampling tt(k) =
%                   (k - 0.5) * gradRasterTime; with oversampling tt(k) =
%                   k * 0.5 * gradRasterTime.
%       .shape_dur  double, seconds, total duration of the shape.
%                   Without oversampling: length(waveform) * gradRasterTime.
%                   With oversampling:  (length(waveform)+1) * 0.5 * gradRasterTime.
%       .first      double, Hz/m, leading boundary value (input or extrapolated)
%       .last       double, Hz/m, trailing boundary value (input or extrapolated)
%
%   ERRORS
%     - 'Slew rate violation (X%)': max(|diff(waveform)|/raster), including
%       the boundary differences against 'first' and 'last', exceeds the
%       active slew limit. Checked before the amplitude check.
%     - 'Gradient amplitude violation (X%)': max(|waveform|) exceeds the
%       active gradient amplitude limit.
%     - 'when oversampling is active the gradient shape vector must contain
%       an odd number of samples': oversampling=true with even-length waveform.
%     - MATLAB:unrecognizedStringChoice: channel is not 'x', 'y', or 'z'.
%
%   NOTES
%     - All gradient values are in Hz/m (Pulseq's internal unit). To convert
%       from mT/m use mr.convert(value, 'mT/m', 'Hz/m', 'gamma', sys.gamma).
%     - The slew-rate check uses the user-supplied (or extrapolated) 'first'
%       and 'last' to compute the boundary slew. Setting 'first' or 'last'
%       far from waveform(1) or waveform(end) can cause spurious slew errors.
%     - Without explicit 'first'/'last' the function emits a deprecation
%       warning on each call; future Pulseq releases will require them.
%     - The 'oversampling' mode is intended for selective-RF gradient shapes
%       where a half-raster sampling reduces aliasing of the encoded waveform.
%     - The check order is slew-then-amplitude, so a wildly out-of-range
%       constant waveform reports as a slew violation rather than an
%       amplitude one.
%
%   EXAMPLE
%     sys = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
%                   'MaxSlew', 170, 'SlewUnit', 'T/m/s');
%     N    = 200;
%     amp  = 0.5e6;          % Hz/m
%     ramp = 80;             % samples in each ramp
%     wf   = [linspace(0, amp, ramp), ...
%             amp * ones(1, N - 2*ramp), ...
%             linspace(amp, 0, ramp)];
%     g = mr.makeArbitraryGrad('x', wf, sys, 'first', 0, 'last', 0);
%     seq = mr.Sequence(sys);
%     seq.addBlock(g);
%
%   SEE ALSO
%     mr.makeTrapezoid, mr.makeExtendedTrapezoid, mr.addGradients,
%     mr.traj2grad, mr.pts2waveform, mr.opts, mr.Sequence/addBlock

persistent parser

if isempty(parser)
    validChannels = {'x','y','z'};
    parser = mr.aux.InputParserCompat;
    parser.FunctionName = 'makeArbitraryGrad';
    parser.addRequired('channel',...
        @(x) any(validatestring(x,validChannels)));
    parser.addRequired('waveform');    
    parser.addOptional('system', [], @isstruct);
    parser.addParamValue('oversampling',false,@islogical);
    parser.addParamValue('maxGrad',0,@isnumeric);
    parser.addParamValue('maxSlew',0,@isnumeric);
    parser.addParamValue('delay',0,@isnumeric);
    parser.addParamValue('first',NaN,@isnumeric);
    parser.addParamValue('last',NaN,@isnumeric);
end
parse(parser,channel,varargin{:});
opt = parser.Results;

if isempty(opt.system)
    system=mr.opts();
else
    system=opt.system;
end

maxSlew=system.maxSlew;
maxGrad=system.maxGrad; % TODO: use this when no duration is supplied
if opt.maxGrad>0
    maxGrad=opt.maxGrad;
end
if opt.maxSlew>0
    maxSlew=opt.maxSlew;
end

g=opt.waveform(:);

if isfinite(opt.first)
    first = opt.first;
else
    warning('it will be compulsory to provide the first point of the gradient shape in the future releases; finding the first by extrapolation for now...');
    if opt.oversampling
        first = 2*g(1)-g(2); % extrapolate by 1 gradient raster
    else
        first = (3*g(1)-g(2))*0.5; % extrapolate by 1/2 gradient of the raster
    end
end

if isfinite(opt.last)
    last = opt.last;
else
    warning('it will be compulsory to provide the last point of the gradient shape in the future releases; finding the last by extrapolation for now...');
    if opt.oversampling
        last = g(end)*2-g(end-1); % extrapolate by 1 gradient raster
    else
        last = (g(end)*3-g(end-1))*0.5; % extrapolate by 1/2 gradient of the raster
    end
end

if opt.oversampling
    slew=[(first-g(1)); (g(2:end)-g(1:end-1)); (last-g(end))]./system.gradRasterTime*2;
else
    slew=[(first-g(1))*2; (g(2:end)-g(1:end-1)); (g(end)-last)*2]./system.gradRasterTime;
end
if ~isempty(slew) && max(abs(slew))>maxSlew
    error('Slew rate violation (%.0f%%)',max(abs(slew))/maxSlew*100);
end
if max(abs(g))>maxGrad
    error('Gradient amplitude violation (%.0f%%)',max(abs(g))/maxGrad*100);
end

grad.type = 'grad';
grad.channel = opt.channel;
grad.waveform = g;
grad.delay = opt.delay;
% true timing and aux shape data
if opt.oversampling
    grad.area=sum(grad.waveform(1:2:end))*system.gradRasterTime; % undo oversamping
    if (mod(length(g),2)~=1)
        error('when oversampling is active the gradient shape vector must contain an odd number of samples');
    end
    grad.tt = (1:length(g))'*0.5*system.gradRasterTime;
    grad.shape_dur = (length(g)+1)*0.5*system.gradRasterTime;    
else
    grad.area=sum(grad.waveform)*system.gradRasterTime; 
    grad.tt = ((1:length(g))'-0.5)*system.gradRasterTime;
    grad.shape_dur = length(g)*system.gradRasterTime;
end
grad.first = first;
grad.last = last;

end

function grad = makeExtendedTrapezoid(channel, varargin)
%makeExtendedTrapezoid Create an extended trapezoid gradient event.
%
%   PURPOSE
%     Build a piecewise-linear gradient event from user-specified
%     (time, amplitude) breakpoints. The breakpoints define the corners
%     of an arbitrary trapezoidal shape (any number of ramps and flats),
%     so this function generalizes mr.makeTrapezoid for cases where the
%     simple rise/flat/fall geometry is not enough (asymmetric ramps,
%     refocusing-and-spoiler combinations, multi-segment readouts, etc.).
%     The returned struct is consumed by mr.Sequence/addBlock to add the
%     gradient to a sequence.
%
%     By default the breakpoint times must already lie on the gradient
%     raster (system.gradRasterTime); this is the "extended trapezoid"
%     mode. Setting 'convert2arbitrary' to true resamples the linear
%     segments onto the raster via mr.pts2waveform and returns an
%     arbitrary-gradient struct, which removes the on-raster requirement
%     for intermediate breakpoints.
%
%   SIGNATURES
%     g = mr.makeExtendedTrapezoid(channel, 'times', t, 'amplitudes', a)
%     g = mr.makeExtendedTrapezoid(channel, system, 'times', t, 'amplitudes', a)
%     g = mr.makeExtendedTrapezoid(channel, 'times', t, 'amplitudes', a, 'system', system)
%     g = mr.makeExtendedTrapezoid(..., 'convert2arbitrary', true)
%     g = mr.makeExtendedTrapezoid(..., 'skip_check', true)
%
%     'times' and 'amplitudes' are required in practice (the parser
%     defaults of 0 trigger the "all times zero" error). Their lengths
%     must match. The first amplitude must be 0 unless this gradient
%     connects to a nonzero gradient in the previous block (in which
%     case set 'skip_check' to true). The last time point must lie on
%     the gradient raster regardless of 'convert2arbitrary'.
%
%   INPUTS
%     channel              [required]    char, gradient axis: 'x', 'y', or 'z'.
%     system               [optional]    struct from mr.opts. May be passed
%                                        positionally OR as 'system', value.
%                                        If omitted, mr.opts() is called for
%                                        defaults.
%     'times'              [name/value]  numeric vector, seconds, breakpoint
%                                        times in strictly ascending order.
%                                        First entry is typically 0; the last
%                                        entry must be on system.gradRasterTime.
%                                        Without 'convert2arbitrary' all entries
%                                        must be on the raster (10 ns rounding
%                                        tolerance). Default 0 (invalid alone).
%     'amplitudes'         [name/value]  numeric vector, Hz/m, gradient
%                                        amplitudes at the corresponding times.
%                                        Same length as 'times'. Default 0.
%     'maxGrad'            [name/value]  numeric, Hz/m, default 0 (use
%                                        system.maxGrad). Only consulted when
%                                        'convert2arbitrary' is true; ignored
%                                        otherwise.
%     'maxSlew'            [name/value]  numeric, Hz/m/s, default 0 (use
%                                        system.maxSlew). Only consulted when
%                                        'convert2arbitrary' is true; ignored
%                                        otherwise.
%     'skip_check'         [name/value]  logical, default false. If true,
%                                        skips the connect-to-previous-block
%                                        check on the first amplitude. Use when
%                                        the gradient resumes a nonzero value
%                                        from the prior block.
%     'convert2arbitrary'  [name/value]  logical, default false. If true,
%                                        resamples the piecewise-linear shape
%                                        onto the gradient raster via
%                                        mr.pts2waveform and dispatches to
%                                        mr.makeArbitraryGrad. Allows
%                                        intermediate breakpoint times to be
%                                        off-raster.
%
%   OUTPUT
%     grad  struct. Field order depends on 'convert2arbitrary'.
%
%     Default (convert2arbitrary = false):
%       .type       char,    always 'grad'
%       .channel    char,    'x' | 'y' | 'z'
%       .waveform   column vector, Hz/m, the input amplitudes
%       .delay      double, seconds, first time point rounded to
%                   gradRasterTime (leading delay before the shape)
%       .tt         column vector, seconds, breakpoint times relative
%                   to .delay (so .tt(1) is 0 or near 0)
%       .shape_dur  double, seconds, last time point minus delay,
%                   rounded to gradRasterTime
%       .area       double, 1/m, trapezoidal integral of the
%                   (tt, waveform) breakpoints (sum of segment areas)
%       .first      double, Hz/m, amplitudes(1)
%       .last       double, Hz/m, amplitudes(end)
%
%     With convert2arbitrary = true (fields delegated to mr.makeArbitraryGrad,
%     so the order changes):
%       .type, .channel, .waveform, .delay, .area, .tt, .shape_dur,
%       .first, .last
%       Here .waveform/.tt/.shape_dur describe the resampled raster
%       waveform from mr.pts2waveform; .area is computed by
%       mr.makeArbitraryGrad. .delay equals times(1) (not rounded here;
%       pts2waveform handles the off-grid first sample). .first/.last
%       are still overwritten with amplitudes(1)/amplitudes(end).
%
%   ERRORS
%     - 'Times and amplitudes must have the same length.': size(times)
%       differs from size(amplitudes).
%     - 'At least one of the given times must be non-zero.': times is
%       all zeros (also fires when 'times' is left at its default 0).
%     - 'Times must be in ascending order and all times must be distinct.':
%       any(diff(times) <= 0).
%     - 'The last time point must be on a gradient raster.': times(end)
%       is not within 10 ns of an integer multiple of gradRasterTime.
%     - 'If first amplitude of a gradient is nonzero, it must connect
%       to previous block!': times(1) > 0 and amplitudes(1) ~= 0 and
%       'skip_check' is false.
%     - 'All time points must be on a gradient raster or
%       "convert2arbitrary" option must be used.': in the default mode,
%       any intermediate time is off-raster (10 ns tolerance).
%     - MATLAB:unrecognizedStringChoice: channel is not 'x', 'y', or 'z'.
%     - When 'convert2arbitrary' is true, mr.makeArbitraryGrad's slew
%       and amplitude checks are inherited (e.g.
%       'Slew rate violation', 'Gradient amplitude violation').
%
%   NOTES
%     - The raster check uses a fixed 10 ns tolerance, independent of
%       gradRasterTime, so very coarse rasters still accept tightly
%       integer-multiple times.
%     - In the default mode the segment-area integrator runs on the
%       relative times .tt, so the reported .area excludes any leading
%       .delay region (which has zero gradient by construction).
%     - .first and .last are always set from amplitudes(1) and
%       amplitudes(end) at the end of the function, even in the
%       convert2arbitrary branch where mr.makeArbitraryGrad would
%       otherwise extrapolate them.
%     - 'maxGrad' and 'maxSlew' overrides only take effect in the
%       convert2arbitrary branch; the default branch performs no slew
%       or amplitude check at all and silently produces shapes that
%       may exceed system limits.
%     - Caches an inputParser in a persistent variable for performance;
%       no other global state.
%
%   EXAMPLE
%     sys = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
%                   'MaxSlew', 170, 'SlewUnit', 'T/m/s');
%     % Asymmetric trapezoid: short rise, long flat, long fall
%     amp = sys.maxGrad * 0.5;
%     t   = [0, 100e-6, 600e-6, 1000e-6];     % all on gradRasterTime
%     a   = [0, amp,    amp,    0];
%     g   = mr.makeExtendedTrapezoid('z', sys, 'times', t, 'amplitudes', a);
%
%     % Same shape but with off-raster intermediate times: needs
%     % convert2arbitrary to resample onto the gradient raster.
%     t2  = [0, 95e-6, 605e-6, 1000e-6];
%     g2  = mr.makeExtendedTrapezoid('z', sys, ...
%               'times', t2, 'amplitudes', a, 'convert2arbitrary', true);
%
%   SEE ALSO
%     mr.makeTrapezoid, mr.makeArbitraryGrad, mr.opts, mr.pts2waveform,
%     mr.calcDuration, mr.Sequence/addBlock
%
%   Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>

persistent parser

if isempty(parser)
    validChannels = {'x', 'y', 'z'};
    parser = mr.aux.InputParserCompat;
    parser.FunctionName = 'makeExtendedTrapezoid';
    parser.addRequired('channel', ...
        @(x) any(validatestring(x, validChannels)));
    parser.addOptional('system',[],@isstruct);
    parser.addParamValue('times', 0, @isnumeric);
    parser.addParamValue('amplitudes', 0, @isnumeric);
    parser.addParamValue('maxGrad', 0, @isnumeric);
    parser.addParamValue('maxSlew', 0, @isnumeric);
    parser.addParamValue('skip_check', false);
    parser.addParamValue('convert2arbitrary', false);
    
end
parse(parser,channel,varargin{:});
opt = parser.Results;

if isempty(opt.system)
    system=mr.opts();
else
    system=opt.system;
end

if any(size(opt.times) ~= size(opt.amplitudes))
    error('Times and amplitudes must have the same length.');
end

if all(opt.times == 0)
    error('At least one of the given times must be non-zero.');
end

if any(diff(opt.times)<=0)
    error('Times must be in ascending order and all times must be distinct.');
end

if abs(round(opt.times(end)/system.gradRasterTime)*system.gradRasterTime-opt.times(end))>1e-8 % 10ns is an acceptable rounding error
    error('The last time point must be on a gradient raster.');
end

%if all(opt.amplitudes == 0)
%    error('At least one of the given amplitudes must be non-zero.');
%end

if opt.skip_check == false && opt.times(1) > 0 && opt.amplitudes(1) ~= 0
    error('If first amplitude of a gradient is nonzero, it must connect to previous block!');
end

maxSlew = system.maxSlew;
maxGrad = system.maxGrad;
if opt.maxGrad > 0
    maxGrad = opt.maxGrad;
end
if opt.maxSlew > 0
    maxSlew = opt.maxSlew;
end

if (opt.convert2arbitrary)
    % represent the extended trapezoid on the regularly sampled time grid
    waveform = mr.pts2waveform(opt.times, opt.amplitudes, system.gradRasterTime);
    grad = mr.makeArbitraryGrad(channel, waveform, system, ...
                                'maxSlew', maxSlew,...
                                'maxGrad', maxGrad,...
                                'delay', opt.times(1));
else
    % keep the original possibly irregular sampling
    if any(abs(round(opt.times/system.gradRasterTime)*system.gradRasterTime-opt.times)>1e-8) % 10ns is an acceptable rounding error
        error('All time points must be on a gradient raster or "convert2arbitrary" option must be used.');
    end
    %
    grad.type = 'grad';
    grad.channel = opt.channel;
    grad.waveform = opt.amplitudes(:);
    grad.delay = round(opt.times(1)/system.gradRasterTime)*system.gradRasterTime;
    grad.tt = opt.times(:) - grad.delay;
    grad.shape_dur = round(grad.tt(end)/system.gradRasterTime)*system.gradRasterTime;
    grad.area=0.5*sum((grad.tt(2:end)-grad.tt(1:end-1)).*(grad.waveform(2:end)+grad.waveform(1:end-1)));
end

% MZ: although makeArbitraryGrad sets the .first and .last for extended 
% trapezoids we can do it better
grad.first=opt.amplitudes(1);
grad.last=opt.amplitudes(end);

end



% figure; plot(waveform)
% waveform = waveform(1:end-2);
% opt.times = round(opt.times/system.gradRasterTime)*system.gradRasterTime; % round onto grid
% times_diff = diff(opt.times);
% amplitudes_diff = diff(opt.amplitudes);
% waveform = [];
% for ii = 1:length(opt.times)-1
%     % SK: there are no new points after the end, therefore we dont need to
%     % handle the overlap situation.
%     if ii == length(opt.times)-1
%         crop = 0;
%     else
%         crop = system.gradRasterTime;
%     end
%     y = amplitudes_diff(ii)/times_diff(ii)*...
%         (0:system.gradRasterTime:(opt.times(ii+1)-opt.times(ii)-crop))...
%         + opt.amplitudes(ii);
%     waveform = [waveform y(1:end)];
% end

function del = makeDelay(delay)
%makeDelay Create a delay event.
%
%   PURPOSE
%     Build a delay event struct of the requested duration. A delay event
%     pads a block to a specified length without producing any RF, gradient,
%     or ADC activity. The returned struct is consumed by
%     mr.Sequence/addBlock, either alone (a pure-delay block) or together
%     with other events to extend the block duration to the delay value.
%
%   SIGNATURES
%     del = mr.makeDelay(delay)        % single positional argument, in seconds
%
%     There is exactly one calling form. The delay is a positional scalar
%     in seconds; there are no name/value parameters and no system argument.
%     When passed to addBlock alongside other events, the resulting block
%     duration is the maximum of the delay and the other events' durations.
%
%   INPUTS
%     delay  [required]  numeric scalar, seconds, must be finite and >= 0.
%                        A value of 0 is accepted (produces a zero-length
%                        delay event, useful as a placeholder).
%
%   OUTPUT
%     del  struct with fields (in order returned by fieldnames):
%       .type   char, always 'delay'
%       .delay  double, seconds, the requested delay duration
%
%   ERRORS
%     - makeDelay:invalidArguments 'Must supply a delay': raised when called
%       with no arguments.
%     - makeDelay:invalidDelay 'Delay (%.2f ms) is invalid': raised when the
%       delay is NaN, Inf, or negative. The message reports the offending
%       value in milliseconds.
%
%   NOTES
%     - The delay is not snapped to any raster here. The caller is
%       responsible for rounding to system.gradRasterTime (or another raster
%       appropriate to the block contents) when timing must align.
%     - In a multi-event block, addBlock takes the maximum of all event
%       durations; a delay shorter than the longest event has no effect on
%       block duration.
%
%   EXAMPLE
%     sys = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
%                   'MaxSlew', 170, 'SlewUnit', 'T/m/s');
%     seq = mr.Sequence(sys);
%     rf  = mr.makeBlockPulse(pi/2, 'Duration', 1e-3, 'system', sys);
%     adc = mr.makeAdc(256, 'Duration', 5e-3, 'system', sys);
%     TE  = 20e-3;
%     delayTE = TE - mr.calcDuration(rf)/2 - mr.calcDuration(adc)/2;
%     delayTE = ceil(delayTE/sys.gradRasterTime)*sys.gradRasterTime;
%     seq.addBlock(rf);
%     seq.addBlock(mr.makeDelay(delayTE));
%     seq.addBlock(adc);
%
%   SEE ALSO
%     mr.Sequence/addBlock, mr.calcDuration, mr.makeSoftDelay

if nargin<1
    error('makeDelay:invalidArguments','Must supply a delay');
end
assert(isfinite(delay) & delay>=0,'makeDelay:invalidDelay',...
    'Delay (%.2f ms) is invalid',delay*1e3);
del.type = 'delay';
del.delay = delay;
end

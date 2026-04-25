function trig = makeDigitalOutputPulse(channel, varargin)
%makeDigitalOutputPulse Create a digital output pulse event a.k.a. trigger.
%
%   PURPOSE
%     Build a digital output trigger event on one of the scanner's hardware
%     output channels. The event holds a single channel high for the
%     requested duration after an optional delay. The returned struct is
%     consumed by mr.Sequence/addBlock and is typically packed alongside
%     gradient/RF/ADC events in the same block (in which case the trigger
%     simply rides along on that block's timeline). Channel meaning is
%     interpreter/scanner specific; the sequence file just records the
%     channel name.
%
%   SIGNATURES
%     trig = mr.makeDigitalOutputPulse(channel)
%     trig = mr.makeDigitalOutputPulse(channel, 'duration', dur)
%     trig = mr.makeDigitalOutputPulse(channel, 'duration', dur, 'delay', d)
%     trig = mr.makeDigitalOutputPulse(channel, ..., 'system', sys)
%
%     channel is the only positional argument and must be one of the three
%     supported names below. All other parameters are name/value pairs.
%
%   INPUTS
%     channel     [required]    char, one of 'osc0', 'osc1', 'ext1'. Selects
%                               which hardware output line the trigger
%                               drives. Any other value raises an error.
%     'delay'     [name/value]  numeric scalar, seconds, default 0. Delay
%                               from the start of the enclosing block until
%                               the trigger goes high.
%     'duration'  [name/value]  numeric scalar, seconds, default 0. Time the
%                               trigger stays high. If the requested value
%                               is <= system.gradRasterTime, it is silently
%                               raised to system.gradRasterTime (so the
%                               default 0 produces a one-raster pulse).
%     'system'    [name/value]  struct from mr.opts, default mr.opts() (the
%                               built-in defaults). Only system.gradRasterTime
%                               is consulted, to enforce the minimum
%                               duration.
%
%   OUTPUT
%     trig  struct with fields (in order returned by fieldnames):
%       .type      char, always 'output'
%       .channel   char, the channel name passed in ('osc0', 'osc1', or 'ext1')
%       .delay     double, seconds, the delay as supplied (not rasterized)
%       .duration  double, seconds, the duration after the minimum-raster bump
%
%   ERRORS
%     - makeDigitalOutputPulse:invalidArguments 'Must supply a channel':
%       raised when called with no arguments.
%     - makeDigitalOutputPulse:invalidChannel 'Channel (%s) is invalid':
%       raised when channel is not one of 'osc0', 'osc1', 'ext1'.
%
%   NOTES
%     - The duration bump is silent: any requested duration <= gradRasterTime
%       (including the default 0) becomes exactly gradRasterTime in the
%       returned struct. To get a longer trigger, pass an explicit duration
%       greater than gradRasterTime.
%     - The delay is not snapped to any raster. The caller is responsible
%       for ensuring the delay aligns with the block's timing constraints.
%     - 'output' triggers are distinct from 'trigger' (input) events created
%       by mr.makeTrigger; output events drive the scanner's external lines,
%       input events wait for an external signal.
%
%   EXAMPLE
%     sys = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
%                   'MaxSlew', 170, 'SlewUnit', 'T/m/s');
%     trig = mr.makeDigitalOutputPulse('osc0', 'duration', 100e-6);
%     seq = mr.Sequence(sys);
%     rf  = mr.makeBlockPulse(pi/2, 'Duration', 1e-3, 'system', sys);
%     seq.addBlock(rf, trig);   % trigger fires alongside the RF pulse
%
%   SEE ALSO
%     mr.Sequence/addBlock, mr.makeTrigger, mr.opts

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeDigitalOutputPulse';
    
    addParamValue(parser, 'delay', 0, @isnumeric);
    addParamValue(parser, 'duration', 0, @isnumeric); % will replace with gradRadterTime below
    addParamValue(parser, 'system', [], @isstruct); 
end

if nargin<1
    error('makeDigitalOutputPulse:invalidArguments','Must supply a channel');
end

parse(parser, varargin{:});
opt = parser.Results;

if isempty(opt.system)
    system=mr.opts();
else
    system=opt.system;
end

channel_num=find(strcmp(channel,{'osc0','osc1','ext1'}));
assert(~isempty(channel_num) && channel_num>0,'makeDigitalOutputPulse:invalidChannel',...
    'Channel (%s) is invalid',channel);
trig.type = 'output';
trig.channel=channel;
trig.delay = opt.delay;
trig.duration = opt.duration;
if (trig.duration<=system.gradRasterTime)
    trig.duration=system.gradRasterTime;
end

end

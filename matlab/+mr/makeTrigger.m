function trig = makeTrigger(channel, varargin)
%makeTrigger Create a trigger halt (input sync) event.
%
%   PURPOSE
%     Build a trigger event struct that causes the scanner sequencer to
%     pause execution and wait for an external input signal on the named
%     channel before resuming. Optionally a pre-sync delay and a post-sync
%     duration can be specified. The returned struct is consumed by
%     mr.Sequence/addBlock to attach the trigger to a sequence block.
%
%     Note: this is an input/halt trigger (sequence waits for a signal),
%     not an output trigger pulse. For output triggers see
%     mr.makeDigitalOutputPulse.
%
%   SIGNATURES
%     trig = mr.makeTrigger(channel)
%     trig = mr.makeTrigger(channel, delay, duration)
%     trig = mr.makeTrigger(channel, delay, duration, system)
%     trig = mr.makeTrigger(channel, 'delay', d, 'duration', dur, 'system', sys)
%
%     The channel argument is required and positional. The delay,
%     duration, and system arguments are positional-or-name/value (the
%     custom mr.aux.InputParserCompat parser accepts both forms even
%     though they are registered as positional). If duration is less than
%     or equal to system.gradRasterTime, it is silently bumped up to
%     system.gradRasterTime.
%
%   INPUTS
%     channel    [required]    char vector, trigger input channel name.
%                              Valid values: 'physio1', 'physio2'
%                              (Siemens-specific). Other values raise
%                              makeTrigger:invalidChannel.
%     delay      [optional]    numeric scalar, seconds. Time from the
%                              start of the block until the sequencer
%                              begins waiting for the trigger signal.
%                              Default 0.
%     duration   [optional]    numeric scalar, seconds. Total length of
%                              the trigger event after the sync. Silently
%                              raised to system.gradRasterTime if smaller.
%                              Default 0 (i.e., bumped to gradRasterTime).
%     system     [optional]    struct, system limits as returned by
%                              mr.opts. Used only to read gradRasterTime
%                              for the silent duration bump described
%                              above. Default mr.opts() (defaults).
%
%   OUTPUT
%     trig  struct with fields (in order returned by fieldnames):
%       .type      char vector, always 'trigger'
%       .channel   char vector, the input channel ('physio1' or 'physio2')
%       .delay     numeric scalar, seconds, pre-sync delay
%       .duration  numeric scalar, seconds, post-sync duration
%                  (>= system.gradRasterTime)
%
%   ERRORS
%     - makeTrigger:invalidArguments: raised if called with no arguments
%       ('Must supply a channel').
%     - makeTrigger:invalidChannel: raised if channel is not 'physio1'
%       or 'physio2' ('Channel (X) is invalid').
%
%   NOTES
%     - duration is silently bumped to system.gradRasterTime when the
%       requested value is <= gradRasterTime. With default mr.opts this
%       bumps any duration <= 1e-5 s up to 1e-5 s.
%     - This event halts the sequencer; the block containing the trigger
%       must reserve enough wall-clock time (delay + duration) for the
%       wait. See mr.calcDuration for block timing.
%     - 'physio1'/'physio2' are Siemens-specific channel labels. Vendors
%       other than Siemens may map these to different physical inputs.
%
%   EXAMPLE
%     % Cardiac-gated block: wait for physio1, then play readout.
%     sys = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
%                   'MaxSlew', 170, 'SlewUnit', 'T/m/s');
%     trig = mr.makeTrigger('physio1', 'duration', 2000e-6);
%     seq = mr.Sequence(sys);
%     seq.addBlock(trig);   % sequencer waits here for the cardiac signal
%
%   SEE ALSO
%     mr.makeDigitalOutputPulse, mr.Sequence/addBlock, mr.opts,
%     mr.calcDuration

persistent parser
if isempty(parser)
    parser = mr.aux.InputParserCompat;
    parser.FunctionName = 'makeTrigger';
    
    addOptional(parser, 'delay', 0, @isnumeric);
    addOptional(parser, 'duration', 0, @isnumeric); % will replace with gradRadterTime below
    addOptional(parser, 'system', [], @isstruct); 
end

if nargin<1
    error('makeTrigger:invalidArguments','Must supply a channel');
end

parse(parser, varargin{:});
opt = parser.Results;

if isempty(opt.system)
    system=mr.opts();
else
    system=opt.system;
end

channel_num=find(strcmp(channel,{'physio1','physio2'}));
assert(~isempty(channel_num) && channel_num>0,'makeTrigger:invalidChannel',...
    'Channel (%s) is invalid',channel);
trig.type = 'trigger';
trig.channel=channel;
trig.delay = opt.delay;
trig.duration = opt.duration;
if (trig.duration<=system.gradRasterTime)
    trig.duration=system.gradRasterTime;
end

end

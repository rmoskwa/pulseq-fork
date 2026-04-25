function sd = makeSoftDelay(varargin)
%makeSoftDelay Create a soft delay extension event.
%
%   PURPOSE
%     Build a soft-delay extension event struct. A soft delay attaches to
%     an otherwise-empty (pure delay) block and instructs the interpreter
%     to rewrite that block's duration at runtime according to:
%         dur = input/factor + offset
%     where 'input' is a user-supplied scalar (e.g. desired TE, TR) given
%     to mr.Sequence/applySoftDelay. Multiple blocks that share the same
%     numID are driven by the same input value, which lets one parameter
%     (TE, TR, ...) control several delay blocks simultaneously. The
%     returned struct is consumed by mr.Sequence/addBlock as the extension
%     payload of a delay block, and later used by applySoftDelay.
%
%   SIGNATURES
%     sd = mr.makeSoftDelay(numID, hint)
%     sd = mr.makeSoftDelay(numID, hint, offset)
%     sd = mr.makeSoftDelay(numID, hint, offset, factor)
%     sd = mr.makeSoftDelay(numID, hint, 'offset', o, 'factor', f)
%
%     numID and hint are positional and required. offset and factor may be
%     supplied positionally or as name/value pairs.
%
%   INPUTS
%     numID     [required]    numeric, non-negative integer identifier shared by all
%                             blocks driven by the same input. The same numID must
%                             carry the same hint string in every block that uses it.
%     hint      [required]    char, non-empty, no whitespace. Human-readable label
%                             (e.g. 'TE', 'TR') written to the .seq file alongside
%                             numID. Blocks sharing a numID must share this hint.
%     offset    [optional]    double, additive offset in the duration equation,
%                             seconds. May be positive or negative. Default 0.
%     'offset'  [name/value]  same as positional offset.
%     factor    [optional]    double, divisor in the duration equation, dimensionless.
%                             May be positive or negative but must be nonzero (zero
%                             is rejected later by mr.Sequence/checkTiming, not here).
%                             Default 1.
%     'factor'  [name/value]  same as positional factor.
%
%   OUTPUT
%     sd  struct with fields (in order returned by fieldnames):
%       .type    char,    always 'softDelay'
%       .num     numeric, copy of numID
%       .hint    char,    copy of hint
%       .offset  double,  seconds
%       .factor  double,  dimensionless
%
%   ERRORS
%     - 'makeSoftDelay: parameter ''hint'' may not contain white space caharacters':
%       hint contains one or more whitespace characters (note: error message
%       contains a typo preserved from the source).
%     - MATLAB:InputParser:ArgumentFailedValidation on 'hint': hint is empty
%       or not a char array.
%     - MATLAB:InputParser:ArgumentFailedValidation on 'numID': numID is not numeric.
%     - MATLAB:InputParser:notEnoughInputs: numID or hint omitted.
%
%   NOTES
%     - The soft delay struct alone does not produce any duration. It must be
%       added to the same block as a pure delay (a block whose only payload is
%       a numeric duration or an mr.makeDelay event), and that block's initial
%       duration sets the default value used if applySoftDelay is never called.
%     - factor == 0 is not caught at construction; it is reported later by
%       mr.Sequence/checkTiming as an invalid soft-delay factor.
%     - All blocks sharing the same numID must yield the same default duration
%       (within ~1e-7 s) when computed as (blockDuration - offset)*factor;
%       otherwise mr.Sequence/checkTiming reports an inconsistency.
%     - Caches an inputParser in a persistent variable for performance;
%       no other global state.
%
%   EXAMPLE
%     % Adjustable TE in a spin-echo block: the empty delay block lasts
%     % delayTE1 seconds by default, but applySoftDelay(seq,'TE',TE_user)
%     % will rewrite it to TE_user/2 + (delayTE1 - TE/2).
%     TE      = 50e-3;
%     delayTE1 = 30e-3;
%     seq = mr.Sequence();
%     seq.addBlock(delayTE1, mr.makeSoftDelay(0, 'TE', ...
%                                             'offset', delayTE1 - TE/2, ...
%                                             'factor', 2));
%
%   SEE ALSO
%     mr.Sequence/addBlock, mr.Sequence/applySoftDelay, mr.makeDelay

persistent parser
if isempty(parser)
    parser = mr.aux.InputParserCompat;
    parser.FunctionName = 'makeSoftDelay';
    
    addRequired(parser, 'numID', @isnumeric);
    addRequired(parser, 'hint', @(x) ischar(x)&&~isempty(x)); 
    addOptional(parser, 'offset', 0, @isnumeric); 
    addOptional(parser, 'factor', 1, @isnumeric); 
end

parse(parser, varargin{:});
opt = parser.Results;

if length(regexp(opt.hint, '(\s+)','split'))>1
    error('makeSoftDelay: parameter ''hint'' may not contain white space characters');
end

if opt.factor==0
    error('makeSoftDelay: parameter ''factor'' must be nonzero');
end

sd.type   = 'softDelay';
sd.num    = opt.numID;
sd.hint   = opt.hint;
sd.offset = opt.offset;
sd.factor = opt.factor;

end

function out = makeLabel(type, label, value)
%makeLabel Create a Pulseq label event.
%
%   PURPOSE
%     Build a label event struct that attaches metadata to the surrounding
%     block when added to a sequence. Labels either set or increment named
%     counters/flags that downstream interpreters (e.g., Siemens, GE) map
%     to MDH fields, parallel-imaging tags, motion-correction hooks, or
%     execution hints. Consumed by mr.Sequence/addBlock; the label itself
%     has no duration and contributes no gradient, RF, or ADC samples.
%
%   SIGNATURES
%     lbl = mr.makeLabel(type, label, value)
%
%     All three arguments are required and positional. There are no
%     name/value forms, no default for any argument, and no system
%     argument. Both 'type' and 'label' are case-sensitive: lowercase is
%     rejected with 'Must supply a valid type' / 'Must supply a valid
%     label'.
%
%   INPUTS
%     type   [required]  char, one of:
%                          'SET'  assign value to the label
%                          'INC'  add value to the running label total
%                                 (increment may be negative)
%     label  [required]  char, the label name. Must be one of the strings
%                        returned by mr.getSupportedLabels(). As of this
%                        version the supported set is:
%                          counters: 'SLC','SEG','REP','AVG','SET','ECO',
%                                    'PHS','LIN','PAR','ACQ','TRID'
%                          flags:    'NAV','REV','SMS','REF','IMA','OFF',
%                                    'NOISE'
%                          control:  'PMC','NOROT','NOPOS','NOSCL','ONCE'
%                        Note: 'SET' appears in both the type domain and
%                        the label domain; the two uses are unrelated.
%     value  [required]  numeric or logical. Integer count (or signed
%                        increment when type=='INC') for counter labels;
%                        true/false (or 0/1) for flag labels; small
%                        non-negative integer for control labels (e.g.,
%                        ONCE takes 0/1/2, TRID takes a positive integer
%                        ID). No range checking is performed here.
%
%   OUTPUT
%     lbl  struct with fields (in order returned by fieldnames):
%       .type   char, internal tag: 'labelset' if type=='SET',
%               'labelinc' if type=='INC'. This is the form
%               mr.Sequence/addBlock dispatches on.
%       .label  char, the label name passed in (uppercase, unchanged).
%       .value  numeric or logical, the value passed in (unchanged;
%               logical true displays as 1).
%
%   ERRORS
%     - 'makeLabel:invalidArguments' / 'Must supply exactly 3 parameters':
%         nargin is not 3.
%     - 'makeLabel:invalidArguments' / 'Must supply a valid label':
%         label is not in mr.getSupportedLabels().
%     - 'makeLabel:invalidArguments' / 'Must supply a valid type':
%         type is not 'SET' or 'INC' (case-sensitive).
%     - 'makeLabel:invalidArguments' / 'Must supply a valid numeric or
%         logical value': value is neither numeric nor logical.
%
%   NOTES
%     - The label event has zero duration. Adding a label-only block
%       between two timed blocks does not advance the sequence clock.
%     - Multiple label events may be added to a single block by passing
%       them all to seq.addBlock (alongside any RF, gradient, ADC, or
%       delay events).
%     - Counters can be reset with type=='SET' and value==0; flags can be
%       cleared with type=='SET' and value==false.
%     - The supported-label list is maintained by mr.getSupportedLabels;
%       check that function for the authoritative current list.
%
%   EXAMPLE
%     % Tag an EPI readout: reset the slice counter, mark the navigator
%     % line, set the k-space center, increment the phase-encode line
%     % counter on each blip, and bump the repetition counter at the end.
%     sys = mr.opts();
%     seq = mr.Sequence(sys);
%     Ny  = 64;
%     % at the start of a slice: reset SLC, set the k-space center
%     seq.addBlock(mr.makeLabel('SET','SLC',0), ...
%                  mr.makeLabel('SET','LIN',round(Ny/2)), ...
%                  mr.makeLabel('SET','NAV',true));
%     % during the readout: increment LIN by 1 each phase blip
%     seq.addBlock(mr.makeDelay(1e-3), mr.makeLabel('INC','LIN',1));
%     % at end of repetition: bump REP
%     seq.addBlock(mr.makeLabel('INC','REP',1));
%
%   SEE ALSO
%     mr.Sequence/addBlock, mr.getSupportedLabels, mr.addCustomLabel


supported_labels=mr.getSupportedLabels();
if nargin~=3
    error('makeLabel:invalidArguments','Must supply exactly 3 parameters');
end
if ~any(ismember(supported_labels,label))
    error('makeLabel:invalidArguments','Must supply a valid label');
end
if ~any(ismember({'SET','INC'},type))
    error('makeLabel:invalidArguments','Must supply a valid type');
end
if ~isnumeric(value) && ~islogical(value)
    error('makeLabel:invalidArguments','Must supply a valid numeric or logical value');
end

switch(type)
    case 'SET'
        out.type = 'labelset';
    case 'INC'
        out.type = 'labelinc';
    otherwise
        error('unknown label type');
end
out.label=label;
out.value=value;
end

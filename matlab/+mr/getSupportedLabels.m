function supported_labels = getSupportedLabels()
%getSupportedLabels Return the list of recognised Pulseq label names.
%
%   PURPOSE
%     Provide the authoritative cell array of label names accepted by
%     mr.makeLabel, mr.Sequence/addBlock, and the .seq read/write paths.
%     The list partitions into three groups: data counters (copied to MDH
%     fields by Siemens interpreters and used by GE), data flags (parallel
%     imaging, navigator, offline, noise tags), and control flags that
%     instruct the interpreter rather than annotating the data. The list
%     is cached in a persistent variable in mr.aux.globalVars so the first
%     call does the assembly and subsequent calls are O(1) lookups.
%
%   SIGNATURES
%     lbls = mr.getSupportedLabels()
%
%     No arguments. No name/value forms. Always returns the full list,
%     including any names previously appended via mr.addCustomLabel during
%     the current MATLAB session.
%
%   INPUTS
%     (none)
%
%   OUTPUT
%     lbls  1xN cell array of char vectors. As shipped, N==23 in this
%           order:
%             counters: 'SLC','SEG','REP','AVG','SET','ECO','PHS','LIN',
%                       'PAR','ACQ','TRID'
%             flags:    'NAV','REV','SMS','REF','IMA','OFF','NOISE'
%             control:  'PMC','NOROT','NOPOS','NOSCL','ONCE'
%           Group meanings:
%             - counters: integer indices for slice, segment, repetition,
%               average, set, echo, phase, line, partition, acquisition,
%               and TR-segment identifier (TRID, used by the GE interpreter
%               to optimise execution).
%             - flags: navigator, reverse readout, simultaneous multi-slice,
%               parallel-imaging reference / imaging, offline (excluded
%               from online recon), noise-adjust scan.
%             - control: motion-correction-aware blocks (PMC); ignore FOV
%               position / rotation / scaling (NOROT, NOPOS, NOSCL); ONCE
%               is a 3-state flag (0 every repeat, 1 first only, 2 last
%               only).
%           Custom names registered via mr.addCustomLabel are appended at
%           the end in the order they were registered.
%
%   SIDE EFFECTS
%     - First call (or first call after mr.aux.globalVars('reset',
%       'SupportedLabels'), 'clear all', or 'clear functions') populates
%       a persistent SupportedLabels variable inside mr.aux.globalVars.
%       Subsequent calls in the same session return the cached value.
%       mr.addCustomLabel mutates this same cache.
%
%   NOTES
%     - Comparison against this list is case-sensitive (mr.makeLabel uses
%       ismember). 'slc' is not 'SLC'.
%     - 'SET' appears here as a label name and is also the type token
%       accepted by mr.makeLabel (alongside 'INC'). The two uses are
%       independent; mr.makeLabel('SET','SET',1) is legal.
%     - This list is the master allow-list for label names. .seq files
%       written by mr.Sequence/write encode labels by their position in
%       this list, so reordering or trimming the built-in entries would
%       break round-tripping.
%
%   EXAMPLE
%     % Print the supported labels and check whether a name is recognised.
%     lbls = mr.getSupportedLabels();
%     fprintf('%d supported labels\n', numel(lbls));
%     disp(lbls);
%     isKnown = any(strcmp(lbls, 'LIN'));   % true
%     isKnown = any(strcmp(lbls, 'lin'));   % false (case-sensitive)
%
%   SEE ALSO
%     mr.makeLabel, mr.addCustomLabel, mr.Sequence/addBlock

supported_labels=mr.aux.globalVars('get','SupportedLabels');
if isempty(supported_labels)
  supported_labels={ ...
    ... % data counters
    'SLC','SEG','REP','AVG','SET','ECO','PHS','LIN','PAR','ACQ', ... % these are copied to the corresponding MDH fields on Siemens
    'TRID', ... # an integer ID of the TR (sequence segment) used by the GE interpreter (and some others) to optimize the execution on the scanner
    ... % data flags
    'NAV','REV','SMS', ... % these are copied to the corresponding MDH fields on Siemens 
    'REF', 'IMA', ... % flags for parallel imaging
    'OFF', ... % Offline flag that labels the data, that should not be used for the online-reconstruction (on Siemens it negates the ONLINE MDH flag)
    'NOISE', ... % flag marking noise adjust scan, for parallel imaging acceleration
    ... % control flags/switches -- they are not affecting the data but rather the sequence itself
        'PMC', ... # for MoCo/PMC Pulseq version to recognize blocks that can be prospectively corrected for motion 
        'NOROT','NOPOS','NOSCL', ... # instruct the interpreter to ignore the position, rotation or scaling of the FOV specified on the UI 
        'ONCE' ... # a 3-state flag that instructs the interpreter to alter the sequence when executing multiple repeats as follows: blocks with ONCE==0 are executed on every repetition; ONCE==1: only the first repetition; ONCE==2: only the last repetition
    };
  mr.aux.globalVars('set','SupportedLabels',supported_labels);
end

end

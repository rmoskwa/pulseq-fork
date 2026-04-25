function addCustomLabel(new_lbl)
%addCustomLabel Register a new custom data label.
%
%   PURPOSE
%     Append a new label name to the package-wide supported-label list so
%     that subsequent calls to mr.makeLabel(type, name, value) accept the
%     name without raising 'Must supply a valid label'. Useful for
%     interpreters or downstream tools that recognize site-specific tags
%     beyond the built-in counters/flags.
%
%   SIGNATURES
%     mr.addCustomLabel(new_lbl)
%
%     One required positional argument. No return value; the effect is on
%     the persistent supported-label list maintained by mr.aux.globalVars.
%
%   INPUTS
%     new_lbl  [required]  char, the name to register. Conventionally
%                          uppercase (e.g., 'MYTAG') to match the built-in
%                          labels, but no case folding or enforcement is
%                          done. The mr.makeLabel lookup is case-sensitive
%                          (ismember), so 'mytag' and 'MYTAG' are distinct
%                          entries.
%
%   OUTPUT
%     (none) addCustomLabel returns nothing; its effect is the side effect
%     described below.
%
%   ERRORS
%     - 'addCustomLabel: new label should be a character sctring':
%         new_lbl is not a char (e.g., numeric, string class, cell,
%         struct). The literal message preserves the source typo
%         'sctring'.
%
%   SIDE EFFECTS
%     - Appends new_lbl to the persistent SupportedLabels variable owned
%       by mr.aux.globalVars. The change persists for the rest of the
%       MATLAB session and affects every subsequent mr.getSupportedLabels
%       and mr.makeLabel call. Reset with
%       mr.aux.globalVars('reset','SupportedLabels'). 'clear all' and
%       'clear functions' also clear it.
%
%   NOTES
%     - If new_lbl is already on the supported list, the function emits
%       the warning 'addCustomLabel: label %s is already known' and still
%       appends a duplicate entry. The duplicate is harmless for
%       mr.makeLabel lookups (ismember match) but visible in
%       mr.getSupportedLabels output.
%     - There is no removeCustomLabel; to undo a registration, reset the
%       entire list as shown above.
%
%   EXAMPLE
%     % Register a site-specific tag and emit it on a block.
%     mr.addCustomLabel('MYTAG');
%     sys = mr.opts();
%     seq = mr.Sequence(sys);
%     seq.addBlock(mr.makeDelay(1e-3), mr.makeLabel('SET','MYTAG',1));
%
%   SEE ALSO
%     mr.makeLabel, mr.getSupportedLabels

if ~ischar(new_lbl)
    error('addCustomLabel: new label should be a character sctring');
end

supported_labels=mr.getSupportedLabels();
if any(ismember(supported_labels,new_lbl))
   warning('addCustomLabel: label %s is already known',new_lbl); 
end

supported_labels{end+1}=new_lbl;
mr.aux.globalVars('set','SupportedLabels',supported_labels);

end

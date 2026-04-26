function [supported_rf_use, short_rf_use] = getSupportedRfUse()
%getSupportedRfUse Return the list of recognised Pulseq RF-use tags.
%
%   PURPOSE
%     Provide the authoritative cell array of RF 'use' tags accepted by
%     the RF constructors (mr.makeBlockPulse, mr.makeSincPulse,
%     mr.makeGaussPulse, mr.makeAdiabaticPulse, mr.makeArbitraryRf,
%     mr.makeSLRpulse) as the 'use' name/value parameter, and used by
%     mr.Sequence/write when serialising the [RF] section of a .seq file.
%     Optionally also returns a parallel char vector of single-letter
%     short codes (first character of each tag) for compact display or
%     legacy formats.
%
%   SIGNATURES
%     uses = mr.getSupportedRfUse()
%     [uses, short] = mr.getSupportedRfUse()
%
%     No arguments. No name/value forms. Always returns the full built-in
%     list; the set is hard-coded and cannot be extended at runtime.
%
%   INPUTS
%     (none)
%
%   OUTPUT
%     uses   1x7 cell array of char vectors, in this fixed order:
%              'excitation'   tip magnetisation away from equilibrium
%              'refocusing'   refocusing pulse (e.g. 180 in spin-echo)
%              'inversion'    inversion pulse (e.g. for IR preparation)
%              'saturation'   saturation pulse (e.g. fat-sat, outer-vol)
%              'preparation'  generic magnetisation-preparation pulse
%              'other'        catch-all for pulses not fitting the above
%              'undefined'    placeholder; the default when 'use' is
%                             omitted on RF constructors
%     short  (optional) 1x7 char vector of single-letter codes formed by
%            taking the first character of each entry of `uses`, in the
%            same order. As shipped: 'erispou'. All seven leading
%            characters are unique. Only computed when nargout > 1.
%
%   NOTES
%     - Comparison against this list is case-sensitive: the RF
%       constructors validate via mr.aux.InputParserCompat with
%       any(validatestring(...)), so 'excitation' is accepted but
%       'Excitation' is not.
%     - The default 'use' for every RF constructor when the parameter is
%       omitted is 'undefined', not 'excitation'. 'undefined' is a
%       legitimate tag, not an error sentinel.
%     - mr.Sequence/write embeds this list verbatim into the [RF] section
%       header of every .seq file (see write.m). Reordering or trimming
%       entries would break .seq round-tripping with existing readers.
%     - This list is hard-coded inside the function body, not cached in
%       mr.aux.globalVars (unlike mr.getSupportedLabels). There is no
%       runtime registration mechanism analogous to mr.addCustomLabel.
%
%   EXAMPLE
%     % Print the supported RF-use tags and check a candidate value.
%     uses = mr.getSupportedRfUse();
%     fprintf('%d supported RF-use tags\n', numel(uses));
%     disp(uses);
%     isKnown = any(strcmp(uses, 'refocusing'));   % true
%     isKnown = any(strcmp(uses, 'Refocusing'));   % false (case-sensitive)
%
%     % Get the short single-letter codes too.
%     [uses, short] = mr.getSupportedRfUse();
%     fprintf('short codes: %s\n', short);          % erispou
%
%   SEE ALSO
%     mr.makeBlockPulse, mr.makeSincPulse, mr.makeGaussPulse,
%     mr.makeAdiabaticPulse, mr.makeArbitraryRf, mr.makeSLRpulse,
%     mr.getSupportedLabels

supported_rf_use={'excitation','refocusing','inversion','saturation','preparation','other','undefined'};
if nargout>1
    short_rf_use=cell2mat(cellfun (@(x) x(1),supported_rf_use,'un',0));
end

end

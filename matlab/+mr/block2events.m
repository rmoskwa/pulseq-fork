function c = block2events(b)
%block2events Convert a block struct or wrapped event list into a flat cell array.
%
%   PURPOSE
%     Internal helper used by routines that accept either a block struct
%     (typically from mr.Sequence/getBlock) or a varargin list of event
%     structs. The function normalizes the input to a single cell array
%     so downstream code can iterate without case analysis. Called from
%     mr.calcDuration, mr.rotate, mr.rotate3D, mr.Sequence/addBlock, and
%     mr.TransformFOV.
%
%   SIGNATURES
%     c = mr.block2events(b)
%
%     b may be:
%       - a block struct (has an .rf field), e.g. from seq.getBlock(n)
%       - a cell array of event structs (the typical varargin form)
%       - a 1x1 cell wrapping any of the above, possibly nested
%       - an event struct or anything else (returned unmodified)
%     Only one block struct may be passed; mixing a block struct with
%     other arguments is rejected by an assert.
%
%   INPUTS
%     b  [required]  block struct, cell array, or any other value. See
%                    SIGNATURES for accepted forms. Passing a function's
%                    own varargin directly is the typical usage.
%
%   OUTPUT
%     c  cell array (or, when the input is neither a block struct nor a
%        cell, the input value unchanged). Contents depend on input:
%        - block struct in: c contains the block's non-empty fields in
%          fieldnames order, namely
%          {blockDuration, rf, gx, gy, gz, adc, ...} with empty fields
%          removed. The leading element is the numeric blockDuration
%          (seconds), not an event struct; callers must tolerate that.
%        - cell array of events in: c is the same cell array, with any
%          enclosing 1x1 cell wrappers stripped iteratively. If the
%          first element of the (post-strip) cell is itself a cell, c
%          is set to that inner cell and any trailing elements of the
%          outer cell are dropped.
%        - any other input in: c = b unchanged.
%
%   ERRORS
%     - 'Only a single block structure can be added': the input contains
%       more than one block struct (a struct array of blocks, or a cell
%       whose first element is a block and whose length is greater than
%       one).
%     - MATLAB:badsubscript 'Index exceeds array bounds': called with
%       an empty cell array; the b{1} access fails.
%
%   NOTES
%     - 1x1 cell wrappers are stripped only while both the outer
%       wrapper has length 1 and its single content is itself a cell.
%       A 1x1 cell wrapping a struct is not stripped; the struct is
%       used as the "first" element to test for a block.
%     - The blockDuration scalar (seconds) is kept in the output cell.
%       Callers that switch on event .type (e.g. mr.calcDuration) must
%       check isnumeric/isstruct before dispatching.
%     - When the input cell's first element is itself a cell, only that
%       inner cell is returned; subsequent outer-cell elements are
%       silently discarded. This branch exists to additionally unwrap
%       inputs of the form {{e1, e2, ...}}.
%
%   EXAMPLE
%     % Direct call: expand a block fetched from a sequence.
%     sys = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
%                   'MaxSlew', 170, 'SlewUnit', 'T/m/s');
%     seq = mr.Sequence(sys);
%     gx  = mr.makeTrapezoid('x', sys, 'Area', 1000, 'Duration', 2e-3);
%     adc = mr.makeAdc(64, sys, 'Duration', 1e-3);
%     seq.addBlock(gx, adc);
%     events = mr.block2events(seq.getBlock(1));
%     % events{1} is the numeric blockDuration; events{2:end} are structs.
%
%   SEE ALSO
%     mr.calcDuration, mr.Sequence/addBlock, mr.Sequence/getBlock,
%     mr.rotate, mr.rotate3D

% strip away 1x1 cell wrapper(s) -- otherwise it conflicts with adding ready-made blocks
while iscell(b) && 1==length(b) && iscell(b{1})
    b=b{1};
end

c = b;    % Assume b is already a cell array of events
if iscell(b)
    first = b{1}; % Use first element to test for block structure
else 
    first = b;
end
if isfield(first, 'rf')
    % Argument is a block structure, copy events to cell array
    % varargin for further processing.
    assert(length(b) == 1, 'Only a single block structure can be added');
%    fields = fieldnames(first)';
%%%%%%%
%     c = {};
%     for f = fields
%         if ~isempty(first.(char(f)))
%             c{end+1} = first.(char(f));
%         end
%     end
%%
%     c = cell(1,length(fields));
%     for i = 1:length(fields)
%         c{i} = first.(char(fields(i)));
%     end
    c=struct2cell(first);
    %c(cellfun(@isempty,c))=[];
    c(cellfun('isempty',c))=[];
    %c=c(~cellfun('isempty',c));
    
    %c(cellfun('isnumeric',c))=[]; % remove blockDuration 

elseif iscell(first)
    c = first;
end

end
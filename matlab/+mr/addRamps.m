function varargout=addRamps(k,varargin)
%addRamps Prepend and append ramp segments to a k-space trajectory.
%
%   PURPOSE
%     Build short k-space segments that ramp from 0 to the start of a
%     supplied trajectory and from its end back to 0 within the gradient
%     and slew-rate limits, then concatenate them onto the trajectory. With
%     'rf', a co-running RF waveform is zero-padded over the same time
%     spans so it stays aligned with the extended trajectory. Used when
%     assembling arbitrary-trajectory excitations or readouts whose start
%     or end gradients are nonzero (e.g. spirals).
%
%   SIGNATURES
%     kout = mr.addRamps(k)                                   % numeric input
%     kout = mr.addRamps(k, system)                           % positional system
%     kout = mr.addRamps(k, 'system', system, ...)            % name/value system
%     [kx,ky,...] = mr.addRamps({kx,ky,...}, ...)             % cell input, one output per channel
%     [kx,ky,...,rf] = mr.addRamps({kx,ky,...}, 'rf', rf, ...)% extend RF shape too
%
%     k may be a numeric matrix [C x N] with C in {1,2,3}, or a cell array
%     {k1,k2,...} of equal-length row vectors. The output shape mirrors the
%     input: numeric in -> single matrix out; cell in -> one row vector per
%     channel in the same order. The ramp lengths at the two ends are set
%     independently by mr.calcRamp and are generally different.
%
%   INPUTS
%     k                  [required]    numeric [C x N] (C in {1,2,3}) or cell
%                                      array of equal-length row vectors.
%                                      Trajectory in k-space units of 1/m,
%                                      sampled on system.gradRasterTime (or
%                                      gradRasterTime/2 when 'gradOversampling'
%                                      is true).
%     system             [optional]    struct from mr.opts. May be passed
%                                      positionally OR as 'system', value. If
%                                      omitted, mr.opts() is called for defaults.
%     'rf'               [name/value]  numeric vector or [], default []. RF
%                                      shape on system.rfRasterTime. When
%                                      supplied, returned as the last output
%                                      with zeros prepended and appended over
%                                      the ramp durations.
%     'maxGrad'          [name/value]  numeric, Hz/m, default 0 (use
%                                      system.maxGrad). Override gradient
%                                      amplitude limit for the ramp calculation.
%     'maxSlew'          [name/value]  numeric, Hz/m/s, default 0 (use
%                                      system.maxSlew). Override slew-rate
%                                      limit for the ramp calculation.
%     'gradOversampling' [name/value]  logical, default false. If true the
%                                      ramps are computed on gradRasterTime/2
%                                      to align with an oversampled trajectory.
%
%   OUTPUT
%     varargout depends on the input form:
%       - numeric k input:  single matrix [C x (Nup + N + Ndown)] in 1/m,
%                           same number of channels as the input.
%       - cell k input:     one row vector per channel, in the order given
%                           in the cell, all of length Nup + N + Ndown.
%       - 'rf' supplied:    a final extra output containing
%                           [zeros(1, Nup*10), rf, zeros(1, Ndown*10)]; the
%                           factor 10 assumes the default raster ratio
%                           rfRasterTime / gradRasterTime = 1us / 10us.
%
%   ERRORS
%     - 'Failed to calculate gradient ramps': mr.calcRamp could not connect
%       0 to k(:,1) (or k(:,end) to 0) within its MaxPoints budget under the
%       active gradient and slew limits. Triggered by very large edge
%       gradients; relax limits, oversample, or shrink the trajectory edge
%       gradients to recover.
%
%   NOTES
%     - The 'rf' zero-pad uses a hardcoded factor of 10 between rfRasterTime
%       and gradRasterTime. With a non-default raster ratio the padded RF
%       length will not match the gradient duration.
%     - Internally the trajectory is padded with zero rows to 3 channels for
%       the ramp calculation, then truncated back to the input channel count.
%       A 1D or 2D trajectory therefore still respects the 3D vector slew
%       limit when the ramps are computed.
%     - The number of ramp samples at start vs end (Nup, Ndown) is generally
%       different when the gradient at the two trajectory edges differs.
%
%   EXAMPLE
%     sys = mr.opts('MaxGrad',32,'GradUnit','mT/m', ...
%                   'MaxSlew',130,'SlewUnit','T/m/s');
%     T = 8e-3; n = 8; kMax = 100;
%     dT = sys.gradRasterTime;
%     tk = 0:dT:T-dT;
%     kx = kMax*(1-tk/T).*cos(2*pi*n*tk/T);
%     ky = kMax*(1-tk/T).*sin(2*pi*n*tk/T);
%     tr = 0:sys.rfRasterTime:T-sys.rfRasterTime;
%     signal = exp(-((1-tr/T)*5).^2);
%     [kx,ky,signal] = mr.addRamps({kx,ky}, 'rf', signal, 'system', sys);
%     gx = mr.makeArbitraryGrad('x', mr.traj2grad(kx), 'first', 0, 'last', 0);
%     gy = mr.makeArbitraryGrad('y', mr.traj2grad(ky), 'first', 0, 'last', 0);
%
%   SEE ALSO
%     mr.calcRamp, mr.makeArbitraryGrad, mr.makeArbitraryRf, mr.traj2grad

persistent parser
if isempty(parser)
    parser = mr.aux.InputParserCompat;
    parser.FunctionName = 'addRamps';
    parser.addRequired('k',@(x)(isnumeric(x)||iscell(x)));
    parser.addOptional('system',[],@isstruct);
    parser.addParamValue('rf',[],@isnumeric);
    parser.addParamValue('maxGrad',0,@isnumeric);
    parser.addParamValue('maxSlew',0,@isnumeric);
    parser.addParamValue('gradOversampling',false,@islogical);
    
end
parse(parser,k,varargin{:});
opt = parser.Results;

if isempty(opt.system)
    system=mr.opts();
else
    system=opt.system;
end

if opt.maxGrad>0
    system.maxGrad=opt.maxGrad;
end
if opt.maxSlew>0
    system.maxSlew=opt.maxSlew;
end

if iscell(opt.k)
    k=cell2mat(opt.k(:));
else
    k=opt.k;
end

nChannels=size(k,1);
k=[k; zeros(3-nChannels,size(k,2))];    % Pad out with zeros if needed

[kUp, ok1]   = mr.calcRamp(zeros(3,2),k(:,1:2),system,'gradOversampling',opt.gradOversampling);
[kDown, ok2] = mr.calcRamp(k(:,end-1:end),zeros(3,2),system,'gradOversampling',opt.gradOversampling);
assert(ok1 & ok2,'Failed to calculate gradient ramps');

kUp = [zeros(3,2), kUp];            % Add start and end points to ramps
kDown = [kDown, zeros(3,1)];

k = [kUp, k, kDown];                % Add ramps to trajectory

if isnumeric(opt.k)
    varargout{1} = k(1:nChannels,:);
else
    for i=1:nChannels
        varargout{i}=k(i,:);
    end
end
if ~isempty(opt.rf)
    varargout{end+1}=[zeros(1,size(kUp,2)*10), opt.rf, zeros(1,size(kDown,2)*10)];
end


end
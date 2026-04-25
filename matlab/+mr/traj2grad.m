function [g sr]=traj2grad(k,varargin)
%traj2grad Convert a k-space trajectory to gradient waveform and slew rate.
%
%   PURPOSE
%     Convert a k-space trajectory k(t) (in 1/m) into the gradient
%     waveform g(t) (in Hz/m) that produces it, and an estimate of the
%     slew rate sr(t) (in Hz/m/s). The trajectory samples are assumed to
%     lie on raster edges; the returned gradient samples are the finite
%     differences of k along time and therefore live on raster cell
%     centers, while the slew rate samples are the time derivative of g
%     and live between gradient points (with the first slew sample
%     stitched to the user-supplied 'first' value, see NOTES).
%     Used to design custom-trajectory readouts (spirals, EPI variants,
%     selective-RF k-space-trajectory pulses) and to check that a
%     candidate trajectory respects the system's gradient-amplitude and
%     slew-rate limits before passing the gradient to mr.makeArbitraryGrad.
%
%   SIGNATURES
%     g          = mr.traj2grad(k)
%     [g, sr]    = mr.traj2grad(k)
%     [g, sr]    = mr.traj2grad(k, 'RasterTime', T)
%     [g, sr]    = mr.traj2grad(k, 'system', sys)
%     [g, sr]    = mr.traj2grad(k, 'first', f, 'firstGradStepHalfRaster', tf, ...
%                               'conservativeSlewEstimate', tf)
%
%     k is the only positional argument; everything else is name/value.
%     The output sr is optional; if requested its size matches g exactly.
%
%   INPUTS
%     k                         [required]    double matrix [nChannel nTime], 1/m. Trajectory
%                                             samples on raster edges. Each row is one logical
%                                             channel (e.g., row 1 = kx, row 2 = ky).
%     'RasterTime'              [name/value]  scalar double, seconds. Time between adjacent k
%                                             samples. Default system.gradRasterTime.
%     'system'                  [name/value]  struct from mr.opts. If omitted, mr.opts() defaults
%                                             are used. Only system.gradRasterTime is read, and
%                                             only when 'RasterTime' is not supplied.
%     'first'                   [name/value]  double column vector [nChannel 1], Hz/m. Gradient
%                                             value just before the first sample, used only to
%                                             compute sr(:,1). Default zeros(nChannel,1).
%     'firstGradStepHalfRaster' [name/value]  logical scalar. If true (default), sr(:,1) is
%                                             multiplied by 2 to account for the half-raster
%                                             offset between the trajectory edge and the gradient
%                                             cell center at the start of the waveform. Set false
%                                             when k samples and gradient samples sit on the same
%                                             raster (e.g., when oversampling is disabled).
%     'conservativeSlewEstimate'[name/value]  logical scalar. If false (default), interior slew
%                                             samples are the average of the two adjacent
%                                             gradient-step slews. If true, they take the larger
%                                             magnitude of the two. Use true when checking against
%                                             system.maxSlew to avoid underestimating peaks on
%                                             non-smooth trajectories.
%
%   OUTPUT
%     g    double matrix [nChannel, nTime-1], Hz/m. Pointwise time
%          derivative of k: g(:,j) = (k(:,j+1) - k(:,j))/RasterTime.
%          Output is one column shorter than k.
%     sr   double matrix [nChannel, nTime-1], Hz/m/s. Estimate of the
%          gradient time derivative. Always returned with the same shape
%          as g; populated as follows:
%            sr(:,1)   = (g(:,1) - first)/RasterTime, doubled if
%                        firstGradStepHalfRaster is true.
%            sr(:,2)   = (g(:,2) - g(:,1))/RasterTime when
%                        firstGradStepHalfRaster is true; otherwise
%                        averaged/max-abs per the flag.
%            sr(:,k)   for k >= 3 (or k >= 2 when firstGradStepHalfRaster
%                        is false): average or max-abs of two adjacent
%                        gradient-step slew estimates, controlled by
%                        conservativeSlewEstimate.
%
%   NOTES
%     - Output is one sample shorter than the input. nTime k samples
%       produce nTime-1 gradient samples, by construction (forward finite
%       difference, no zero-padding). Downstream callers that need a
%       waveform of length nTime should pad or extend manually.
%     - 'first' affects only sr, not g. Set it to the gradient value
%       active just before the first k sample so the initial slew
%       estimate is meaningful; leave at zero if the trajectory really
%       does start from rest.
%     - The half-raster doubling (firstGradStepHalfRaster=true) reflects
%       Pulseq's standard convention that an arbitrary gradient's first
%       sample sits half a raster cell after the trajectory's first
%       point. When the calling code already places gradient and
%       trajectory samples on the same grid (e.g., when oversampling is
%       disabled in writeSpiral.m: 'firstGradStepHalfRaster',
%       ~gradOversampling), set this flag to false.
%     - conservativeSlewEstimate=true is recommended whenever the result
%       will be compared against system.maxSlew. The averaged estimator
%       hides slew-rate spikes on non-smooth trajectories.
%     - The function computes finite differences only. No raster
%       alignment, smoothing, or amplitude-limit enforcement is
%       performed; that is the caller's responsibility (see
%       writeSpiral.m for a complete trajectory-design loop).
%
%   EXAMPLE
%     % Check whether a candidate spiral-out readout trajectory respects
%     % the system's gradient-amplitude and slew-rate limits before
%     % handing the gradient to mr.makeArbitraryGrad. This mirrors the
%     % feasibility-check pattern used in demoSeq/writeSpiral.m.
%     sys = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
%                   'MaxSlew', 170, 'SlewUnit', 'T/m/s');
%     T       = 10e-3;                         % readout duration, seconds
%     nTurns  = 4;
%     kMax    = 100;                           % 1/m
%     dt      = sys.gradRasterTime;
%     t       = 0:dt:T-dt;
%     kx      = kMax*(t/T).*cos(2*pi*nTurns*t/T);
%     ky      = kMax*(t/T).*sin(2*pi*nTurns*t/T);
%     k       = [kx; ky];
%     [g, sr] = mr.traj2grad(k, 'system', sys, 'conservativeSlewEstimate', true);
%     fprintf('peak |g|  = %.1f%% of maxGrad\n',  100*max(abs(g(:)))/sys.maxGrad);
%     fprintf('peak |sr| = %.1f%% of maxSlew\n', 100*max(abs(sr(:)))/sys.maxSlew);
%
%   SEE ALSO
%     mr.makeArbitraryGrad, mr.opts, mr.calcDuration, mr.Sequence/addBlock

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'traj2grad';
    parser.addRequired('k',@isnumeric);
    parser.addParamValue('first',[],@isnumeric);
    parser.addParamValue('firstGradStepHalfRaster',true,@islogical);
    parser.addParamValue('conservativeSlewEstimate',false,@islogical);    
    parser.addParamValue('system',[],@isstruct);
    parser.addParamValue('RasterTime',[],@isnumeric);
end
parse(parser,k,varargin{:});
opt = parser.Results;
if isempty(opt.system)
    opt.system=mr.opts();
end
if isempty(opt.RasterTime)
    opt.RasterTime=opt.system.gradRasterTime;
end
if isempty(opt.first)
    opt.first=zeros(size(k,1),1); % QC: if the first gradient point is not given, set it to zero. 2025.01.03
end

% Compute finite difference for gradients in Hz/m
%g=([k(:,2:end)-k(:,1:end-1) zeros(size(k,1),1)])/opt.RasterTime; % MZ: with zero-padding
g=[(k(:,2:end)-k(:,1:end-1))/opt.RasterTime]; % MZ: no zero-padding!

% Compute the slew rate (time derivative of the gradient)
sr0=(g-[opt.first g(:,1:end-1)])/opt.RasterTime;
if opt.firstGradStepHalfRaster
    sr0(:,1)=sr0(:,1)*2; % account for the half-step in the beginning of the shape
end

% now we think how to post-process the results 
% gradient is now sampled between the k-points (on raster cell centers)
% whilst the slew rate is between the gradient points, except of the first
% point, which relies on the opt.first value (and may be a bit off anyway,
% but this is the best estimate that we have)
sr=zeros(size(sr0));
sr(:,1)=sr0(:,1);
if (opt.conservativeSlewEstimate)
    if opt.firstGradStepHalfRaster
        sr(:,2)=sr0(:,2);
        sr(:,3:end)=max_abs(sr0(:,2:end-1),sr0(:,3:end));
    else
        sr(:,2:end)=max_abs(sr0(:,1:end-1),sr0(:,2:end));
    end
else
    if opt.firstGradStepHalfRaster
        sr(:,2)=sr0(:,2);
        sr(:,3:end)=0.5*(sr0(:,2:end-1)+sr0(:,3:end));
    else
        sr(:,2:end)=0.5*(sr0(:,1:end-1)+sr0(:,2:end));
    end
end

end

function out=max_abs(in1, in2)
    if size(in2)~=size(in2)
        error('arrays of incompatible sizes');
    end
    abs1gtoe=abs(in1)>=abs(in2);
    out=in1.*abs1gtoe+in2.*(~abs1gtoe);
end


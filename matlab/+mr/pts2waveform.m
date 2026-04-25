function waveform = pts2waveform(times, amplitudes, gradRasterTime)
%pts2waveform Resample piecewise-linear control points onto a uniform raster.
%
%   PURPOSE
%     Convert a piecewise-linear waveform defined by control points
%     (times, amplitudes) into a uniformly sampled vector on a raster of
%     spacing gradRasterTime. Each output sample is the linear interpolant
%     of the control points evaluated at the midpoint of one raster cell.
%     Used internally by mr.makeExtendedTrapezoid (when 'convert2arbitrary'
%     is set) and by mr.addGradients to bring trapezoids and arbitrary
%     gradients onto a common raster before summation. The returned vector
%     is suitable for handing to mr.makeArbitraryGrad.
%
%   SIGNATURES
%     waveform = mr.pts2waveform(times, amplitudes, gradRasterTime)
%
%     All three arguments are required and positional. There are no
%     name/value options.
%
%   INPUTS
%     times           [required]  double vector, seconds. Control-point times. Must
%                                 be monotonically increasing and the same length as
%                                 amplitudes. Need not lie on the raster: endpoints
%                                 are rounded to the nearest gradRasterTime when the
%                                 output grid is built.
%     amplitudes      [required]  double vector, same units as the desired waveform
%                                 (typically Hz/m for gradients). Same length as times.
%     gradRasterTime  [required]  scalar double, seconds. Raster spacing (typically
%                                 system.gradRasterTime, e.g., 10e-6).
%
%   OUTPUT
%     waveform  Nx1 column double, same units as amplitudes. Length is
%               N = round(max(times)/gradRasterTime) - round(min(times)/gradRasterTime).
%               Sample k (1-based) is the linear interpolant of
%               (times, amplitudes) at time
%                 (round(min(times)/gradRasterTime) + k - 0.5)*gradRasterTime,
%               i.e. the midpoint of raster cell k. Always returned as a
%               column vector regardless of the orientation of times and
%               amplitudes.
%
%   NOTES
%     - Endpoints are rounded to the raster, not snapped down or up
%       deterministically. Off-raster control points between the first and
%       last are honored exactly by the linear interpolation; only the
%       extent of the output grid is rounded.
%     - The final raster point is dropped (grd(1:end-1)) so the output has
%       N samples, not N+1. This matches Pulseq's convention that an
%       arbitrary-gradient waveform of N samples occupies N raster cells.
%     - Sampling at cell midpoints (offset gradRasterTime/2 from the grid
%       start) means the output values are not equal to the control-point
%       amplitudes even when a control point lies on the raster. For a
%       linear ramp from 0 to A over the first cell, sample 1 is A/2.
%     - times and amplitudes must have the same length and times must be
%       strictly increasing. Both conditions are enforced by the
%       underlying interp1 call, which throws if violated.
%
%   EXAMPLE
%     % Resample an extended-trapezoid control-point spec onto the
%     % gradient raster, then wrap it as an arbitrary gradient.
%     sys = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
%                   'MaxSlew', 170, 'SlewUnit', 'T/m/s');
%     amp        = mr.convert(20, 'mT/m', 'Hz/m');
%     times      = [0, 200e-6, 400e-6, 600e-6];   % seconds
%     amplitudes = [0, amp, amp, 0];               % Hz/m
%     waveform   = mr.pts2waveform(times, amplitudes, sys.gradRasterTime);
%     g = mr.makeArbitraryGrad('x', waveform, sys);
%
%   SEE ALSO
%     mr.makeExtendedTrapezoid, mr.makeArbitraryGrad, mr.addGradients,
%     mr.opts, interp1

grd = (round(min(times)/gradRasterTime):round(max(times)/gradRasterTime))*gradRasterTime; % the previous code was clipping the gradient now and then...
grd = grd(1:end-1).';
waveform = interp1(times, amplitudes, grd + gradRasterTime/2);

% % times = ceil(times/gradRasterTime)*gradRasterTime; % round onto grid
% times = ceil(times/gradRasterTime); % round onto grid. SK: Dont multiply by 
%                                     % gradRasterTime here: This will 
%                                     % introduce numerical inaccuracies.
% times_diff = diff(times);
% amplitudes_diff = diff(amplitudes);
% waveform = [];
% for ii = 1:length(times)-1
%     % SK: there are no new points after the end, therefore we dont need to
%     % handle the overlap situation.
%     if ii == length(times)-1
%         crop = 0;
%     else
%         crop = 1;
%     end
%     y = amplitudes_diff(ii)/times_diff(ii)*...
%         (0:1:(times(ii+1)-times(ii)-crop))...
%         + amplitudes(ii);
%     waveform = [waveform y(1:end)];
% end
end

function [grad] = scaleGrad(grad, scale, system)
%scaleGrad Scale a gradient event by a scalar.
%
%   PURPOSE
%     Multiplies a gradient event's amplitude (or waveform samples) by a
%     scalar factor and returns a new gradient struct of the same type.
%     Common uses are reversing polarity (scale = -1) and stepping through
%     phase-encoding amplitudes inside a sequence loop. The returned
%     struct can be passed directly to mr.Sequence/addBlock.
%
%     If the optional system limits struct is provided, the scaled
%     gradient is checked against system.maxGrad and system.maxSlew before
%     being returned; out-of-spec results raise an error rather than being
%     silently truncated.
%
%   SIGNATURES
%     g = mr.scaleGrad(grad, scale)            % scale only, no limits check
%     g = mr.scaleGrad(grad, scale, system)    % scale and verify against limits
%
%     The third argument is strictly positional. The 'system', sys
%     name/value form is NOT accepted by this function (unlike most other
%     mr.make* functions).
%
%   INPUTS
%     grad    [required]    gradient event struct from mr.makeTrapezoid,
%                           mr.makeArbitraryGrad, or mr.makeExtendedTrapezoid.
%     scale   [required]    double, dimensionless multiplier. Use -1 for
%                           polarity reversal; any real value is allowed
%                           (limits are checked only if system is provided).
%     system  [optional]    struct from mr.opts(). If omitted, no limits
%                           check is performed. If provided, the scaled
%                           gradient must satisfy system.maxGrad and
%                           system.maxSlew.
%
%   OUTPUT
%     grad struct. Same shape and field order as the input. Field changes
%     depend on the input type:
%
%       Trapezoid input (grad.type == 'trap'):
%         .amplitude   Hz/m, scaled by scale
%         .area        1/m,  scaled by scale
%         .flatArea    1/m,  scaled by scale
%         (.type, .channel, .riseTime, .flatTime, .fallTime, .delay,
%          .first, .last are unchanged)
%
%       Arbitrary or extended trapezoid input (grad.type == 'grad'):
%         .waveform    Hz/m, every sample multiplied by scale
%         .first       Hz/m, scaled by scale
%         .last        Hz/m, scaled by scale
%         (.type, .channel, .delay, .tt, .shape_dur, and importantly
%          .area are NOT updated; see NOTES)
%
%   ERRORS
%     - 'attempting to scale readily registered object! ...': raised when
%       the input grad has an 'id' field, meaning it has already been
%       registered with a Sequence's event library. Either scale the
%       gradient before registering it with addBlock, or deregister first
%       by calling rmfield(grad, 'id').
%     - 'mr.scaleGrad: maximum amplitude exceeded (X %)': scaled gradient
%       amplitude exceeds system.maxGrad. Only raised when system is provided.
%     - 'mr.scaleGrad: maximum slew rate exceeded (X %)': scaled gradient
%       slew rate exceeds system.maxSlew. Only raised when system is provided.
%       For trap, slew is computed against min(riseTime, fallTime); for
%       arbitrary/extended, it is the largest |diff(waveform)/diff(tt)|.
%
%   NOTES
%     - For arbitrary and extended trapezoid inputs, the .area field is
%       NOT recomputed. Only .waveform, .first, and .last are scaled. If
%       you rely on .area downstream after scaling such a gradient,
%       recompute it yourself from the scaled waveform.
%     - Limits checking is opt-in: omitting system skips both the
%       amplitude and slew-rate checks entirely, allowing unbounded
%       scaling. Pass system whenever the scaled gradient is intended for
%       actual playout.
%     - For trap inputs, the slew check uses min(riseTime, fallTime) to
%       find the steepest ramp; both ramps inherit the scaled amplitude.
%
%   EXAMPLE
%     sys = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
%                   'MaxSlew', 170, 'SlewUnit', 'T/m/s');
%     Nx = 64;  fov = 256e-3;  deltak = 1/fov;
%     gx    = mr.makeTrapezoid('x', sys, 'FlatArea', Nx*deltak, ...
%                              'FlatTime', 3.2e-3);
%     gyPre = mr.makeTrapezoid('y', sys, 'Area', deltak*Nx/2, ...
%                              'Duration', 1e-3);
%     % Reverse polarity of the readout for an EPI-style flyback:
%     gx_rev = mr.scaleGrad(gx, -1);
%     % Step phase-encoding amplitude across a loop:
%     peScales = ((0:Nx-1) - Nx/2) / (Nx/2);
%     for i = 1:Nx
%         gy_i = mr.scaleGrad(gyPre, peScales(i), sys);
%         % seq.addBlock(gx, gy_i, ...);
%     end
%
%   SEE ALSO
%     mr.makeTrapezoid, mr.makeArbitraryGrad, mr.makeExtendedTrapezoid,
%     mr.addGradients, mr.rotate, mr.Sequence/addBlock

    if isfield(grad,'id')
        error('attempting to scale readily registered object! please register objects after calling this function or deregister the argument by calling rmfield(...,''id'')');
    end

    if strcmp(grad.type,'trap')
        grad.amplitude=grad.amplitude*scale;
        grad.area=grad.area*scale;
        grad.flatArea=grad.flatArea*scale;
        if nargin>2
            if system.maxGrad<abs(grad.amplitude)
                error("mr.scaleGrad: maximum amplitude exceeded (%g %%)",100*abs(grad.amplitude)/system.maxGrad);
            end
            if abs(grad.amplitude)>eps && system.maxSlew<abs(grad.amplitude)/min(grad.riseTime,grad.fallTime)
                error("mr.scaleGrad: maximum slew rate exceeded (%g %%)",100*abs(grad.amplitude)/min(grad.riseTime,grad.fallTime)/system.maxSlew);
            end
        end
    else
        grad.waveform=grad.waveform*scale;
        grad.first=grad.first*scale;
        grad.last=grad.last*scale;
        if nargin>2
            if system.maxGrad<max(abs(grad.waveform))
                error("mr.scaleGrad: maximum amplitude exceeded (%g %%)",100*max(abs(grad.waveform))/system.maxGrad);
            end
            if max(abs(grad.waveform))>eps
                grad_max_abs_slew=max(abs(diff(grad.waveform)./diff(grad.tt)));
                if system.maxSlew<grad_max_abs_slew
                    error("mr.scaleGrad: maximum slew rate exceeded (%g %%)",100*grad_max_abs_slew/system.maxSlew);
                end
            end
        end
    end
end


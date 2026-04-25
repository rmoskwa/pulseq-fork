function [adcSegments,adcSamplesPerSegment] = calcAdcSeg(numSamples,dwell,system,mode)
%calcAdcSeg Split a long ADC into raster-aligned segments.
%
%   PURPOSE
%     Some scanners (notably Siemens) cannot execute a single ADC event
%     longer than a hardware sample limit (e.g., 8192 samples on Siemens).
%     This helper computes how to split a desired ADC into N equal
%     segments such that:
%       - each segment has at most system.adcSamplesLimit samples,
%       - the per-segment sample count is divisible by
%         system.adcSamplesDivisor (Siemens requirement),
%       - each segment duration is aligned to system.gradRasterTime.
%     The returned product adcSegments*adcSamplesPerSegment is the
%     actual ADC sample count to pass to mr.makeAdc; it can differ from
%     the requested numSamples (smaller for mode 'shorten', larger for
%     'lengthen').
%
%   SIGNATURES
%     [adcSegments, adcSamplesPerSegment] = mr.calcAdcSeg(numSamples, dwell, system)
%     [adcSegments, adcSamplesPerSegment] = mr.calcAdcSeg(numSamples, dwell, system, mode)
%
%     mode is positional only (no name/value form). Default 'shorten'.
%     Valid modes: 'shorten' | 'lengthen'.
%
%   INPUTS
%     numSamples  [required]  double, desired total ADC samples, integer >0
%     dwell       [required]  double, ADC dwell time, seconds; must be an integer
%                             multiple of system.adcRasterTime
%     system      [required]  struct from mr.opts; reads .adcRasterTime,
%                             .gradRasterTime, .adcSamplesLimit, .adcSamplesDivisor
%     mode        [optional]  char, 'shorten' rounds the total sample count down;
%                             'lengthen' rounds up. Default 'shorten'.
%
%   OUTPUT
%     adcSegments           double, integer in [1, 128]. Number of equal
%                           ADC segments. 1 if system.adcSamplesLimit<=0.
%     adcSamplesPerSegment  double, integer >0. Samples per segment;
%                           pass adcSegments*adcSamplesPerSegment as the
%                           numSamples argument to mr.makeAdc.
%
%   ERRORS
%     - 'In mr.calcAdcSeg(...,mode) mode should be either ''shorten'' or
%       ''lengthen''': mode is not one of the two accepted strings.
%     Assertion failures (no message ID, MATLAB:assertion:failed):
%     - system.gradRasterTime is not an integer multiple of
%       system.adcRasterTime (tolerance 1e-9 s).
%     - dwell is not an integer multiple of system.adcRasterTime
%       (tolerance 1e-9 s).
%     - No segmentation found that satisfies
%       adcSamplesPerSegment <= adcSamplesLimit and adcSegments <= 128
%       within the search range. Adjust dwell or try the other mode.
%
%   NOTES
%     - Early exit: if system.adcSamplesLimit <= 0 the function returns
%       adcSegments = 1, adcSamplesPerSegment = numSamples without any
%       raster checks. This is the default in mr.opts (limit defaults to 0).
%     - Hard cap: adcSegments is constrained to <= 128.
%     - The returned total sample count adcSegments*adcSamplesPerSegment
%       can differ from numSamples. For spiral/EPI readouts, recompute
%       trajectory or readout duration from the returned product.
%
%   EXAMPLE
%     sys = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
%                   'MaxSlew', 170, 'SlewUnit', 'T/m/s', ...
%                   'adcSamplesLimit', 8192);
%     % Spiral-style: pick a desired sample count, round dwell to ADC raster
%     adcDwell  = round(2e-6 / sys.adcRasterTime) * sys.adcRasterTime;
%     numWanted = 12345;
%     [nSeg, sPerSeg] = mr.calcAdcSeg(numWanted, adcDwell, sys);
%     adcSamples = nSeg * sPerSeg;
%     adc = mr.makeAdc(adcSamples, sys, 'Dwell', adcDwell);
%
%   SEE ALSO
%     mr.opts, mr.makeAdc, mr.calcDuration, mr.Sequence/addBlock

if system.adcSamplesLimit<=0
    adcSamplesPerSegment=numSamples;
    adcSegments=1;
    return;
end

if ~exist('mode', 'var')
    mode='shorten';
end

if ~strcmp(mode,'shorten') && ~strcmp(mode,'lengthen')
    error('In mr.calcAdcSeg(...,mode) mode should be either ''shorten'' or ''lengthen''');
end

t_eps=1e-9; % TODO: shift it to the system parameters???

iGR=round(system.gradRasterTime/system.adcRasterTime);
assert(abs(system.gradRasterTime/system.adcRasterTime-iGR)<t_eps);

iDwell=round(dwell/system.adcRasterTime);
assert(abs(dwell/system.adcRasterTime-iDwell)<t_eps);

iCommon=lcm(iGR,iDwell); % least common multiplier
samplesStep=iCommon/iDwell;

% Siemens-specific: number of samples should be divisible by system.adcSamplesDivisor
gcd_adcDiv=gcd(samplesStep,system.adcSamplesDivisor);
if gcd_adcDiv~=system.adcSamplesDivisor
    samplesStep=samplesStep*system.adcSamplesDivisor/gcd_adcDiv;
end

if strcmp(mode,'shorten')
    numSamplesStepped=floor(numSamples/samplesStep);
else
    numSamplesStepped=ceil(numSamples/samplesStep);
end

while numSamplesStepped>0 && numSamplesStepped<2*numSamples/samplesStep
    adcSegmentFactors=factor(numSamplesStepped);
    adcSegments=1;        
    if(length(adcSegmentFactors)>1) 
         % we try all permutations and pick the smallest number of segments
        adcSegmentFactorsPerm=perms(adcSegmentFactors);
        adcSegmentFactorsPermProd=cumprod(adcSegmentFactorsPerm');
        adcSegmentCandidates=unique(adcSegmentFactorsPermProd(:)); % this sorts the sequence
        for i=1:(length(adcSegmentCandidates)-1)
            adcSegments=adcSegmentCandidates(i);
            adcSamplesPerSegment=numSamplesStepped*samplesStep/adcSegments;
            if (adcSamplesPerSegment<=system.adcSamplesLimit && adcSegments<=128)
                break
            end
        end
    else
        adcSamplesPerSegment=numSamplesStepped*samplesStep;
    end
    if (adcSamplesPerSegment<=system.adcSamplesLimit && adcSegments<=128)
        break
    end
    if strcmp(mode,'shorten')
        numSamplesStepped=numSamplesStepped-1; % try again with a smaller number of samples
    else
        numSamplesStepped=numSamplesStepped+1; % try again with a greater number of samples
    end
        
end
assert(numSamplesStepped>0); % we could not find a suitable segmentation...
assert(adcSamplesPerSegment>0); 
assert(adcSegments<=128);
end


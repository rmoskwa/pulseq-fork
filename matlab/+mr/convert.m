function out=convert(in,varargin)
%convert Convert numeric values between MR unit conventions.
%
%   PURPOSE
%     Converts gradient amplitudes, slew rates, or RF (B1) amplitudes
%     between physical units (mT/m, T/m/s, uT, ...) and the Hz-based
%     internal units used by Pulseq (Hz/m, Hz/m/s, Hz). Conversions
%     between Tesla-based and Hz-based units use the gyromagnetic
%     ratio (default 42.576 MHz/T for 1H); override via 'gamma' for
%     other nuclei.
%
%     mr.opts calls mr.convert internally on its 'maxGrad', 'maxSlew',
%     and 'maxB1' inputs (controlled by 'gradUnit', 'slewUnit',
%     'b1Unit'), so values passed to mr.opts do not need to be
%     pre-converted. Use mr.convert directly when displaying Pulseq
%     internal values in physical units (e.g., in a report) or when
%     passing physical-unit values into a function that expects
%     Hz-based input.
%
%   SIGNATURES
%     out = mr.convert(in, fromUnit)             % toUnit defaults to category standard
%     out = mr.convert(in, fromUnit, toUnit)     % explicit target unit
%     out = mr.convert(..., 'gamma', G)          % override gyromagnetic ratio (Hz/T)
%
%     If toUnit is omitted (or passed as []), the result is in the
%     Hz-based "standard" unit for the same category as fromUnit:
%       B1 amplitude  ->  'Hz'
%       gradient      ->  'Hz/m'
%       slew rate     ->  'Hz/m/s'
%     'in' may be a scalar, vector, matrix, or N-D array; the output
%     has the same shape.
%
%   INPUTS
%     in        [required]    numeric, any shape, value(s) to convert
%     fromUnit  [required]    char, one of the unit strings below
%     toUnit    [optional]    char, one of the unit strings below; default =
%                             category standard
%     gamma     [name/value]  numeric, gyromagnetic ratio in Hz/T; default
%                             42.576e6 (1H)
%
%     Valid unit strings, grouped by category:
%       B1 amplitude:   'Hz', 'T', 'mT', 'uT'
%       gradient:       'Hz/m', 'mT/m', 'rad/ms/mm'
%       slew rate:      'Hz/m/s', 'mT/m/ms', 'T/m/s', 'rad/ms/mm/ms'
%
%   OUTPUT
%     out  numeric, same shape as in, expressed in toUnit
%
%   ERRORS
%     - 'MATLAB:unrecognizedStringChoice': fromUnit or toUnit is not
%       one of the 11 valid unit strings (matched via validatestring).
%     - 'MATLAB:InputParser:ArgumentFailedValidation': in is not numeric.
%
%   NOTES
%     - No cross-category check is performed. Passing fromUnit and
%       toUnit from different categories (e.g., 'mT/m' to 'Hz') does
%       not error, but produces a numerically meaningless result
%       because the conversion routes through a single Hz-based
%       standard value. Always pair fromUnit and toUnit within the
%       same category.
%     - The default gamma (42.576e6 Hz/T) is the 1H gyromagnetic ratio.
%       For other nuclei, pass the appropriate gamma via name/value; 
%       otherwise gradient and B1 conversions involving Tesla-based 
%       units will be wrong by the ratio of gyromagnetic ratios.
%     - 'rad/ms/mm' (gradient) and 'rad/ms/mm/ms' (slew) do not depend
%       on gamma; they are pure angular-frequency conversions
%       (factor of 2*pi and a power of ten relative to Hz-based units).
%
%   EXAMPLE
%     % Convert physical-unit limits to Pulseq's internal Hz-based units.
%     gradHz = mr.convert(40,  'mT/m',  'Hz/m');    % gradient amplitude
%     slewHz = mr.convert(170, 'T/m/s', 'Hz/m/s');  % slew rate
%     b1Hz   = mr.convert(20,  'uT',    'Hz');      % RF amplitude
%
%     % Display Pulseq-internal values back in physical units (e.g., for a report).
%     fprintf('Max gradient: %.2f mT/m\n',  mr.convert(gradHz, 'Hz/m',   'mT/m'));
%     fprintf('Max slew:     %.2f T/m/s\n', mr.convert(slewHz, 'Hz/m/s', 'T/m/s'));
%
%     % Non-proton nucleus: 13C (gyromagnetic ratio 10.708 MHz/T).
%     gradHz_13C = mr.convert(40, 'mT/m', 'Hz/m', 'gamma', 10.708e6);
%
%   SEE ALSO
%     mr.opts

persistent parser
validB1Units={'Hz','T','mT','uT'}; % todo: gauss?
validGradUnits={'Hz/m','mT/m','rad/ms/mm'};
validSlewUnits={'Hz/m/s','mT/m/ms','T/m/s','rad/ms/mm/ms'};
validUnits=cat(2,validB1Units,validGradUnits,validSlewUnits);
if isempty(parser)
    parser = mr.aux.InputParserCompat;
    parser.FunctionName = 'convert';
    parser.addRequired('in',@isnumeric);
    parser.addRequired('fromUnit',...
        @(x) any(validatestring(x,validUnits)));
    parser.addOptional('toUnit',[],...
        @(x) any(validatestring(x,validUnits)));
    parser.addParamValue('gamma',42.576e6,@isnumeric); % Hz/T
end
parse(parser,in,varargin{:});
opt = parser.Results;

% Set default output unit if not given
if isempty(opt.toUnit)
    if ismember(opt.fromUnit,validGradUnits)
        opt.toUnit = validGradUnits{1};
    elseif ismember(opt.fromUnit,validSlewUnits)
        opt.toUnit = validSlewUnits{1};
    elseif ismember(opt.fromUnit,validB1Units)
        opt.toUnit = validB1Units{1};
    end
end

% Verify fromUnit and toUnit are in the same category.
if ismember(opt.fromUnit,validB1Units),       fromCat = 'B1';
elseif ismember(opt.fromUnit,validGradUnits), fromCat = 'gradient';
elseif ismember(opt.fromUnit,validSlewUnits), fromCat = 'slew rate';
end

if ismember(opt.toUnit,validB1Units),         toCat = 'B1';
elseif ismember(opt.toUnit,validGradUnits),   toCat = 'gradient';
elseif ismember(opt.toUnit,validSlewUnits),   toCat = 'slew rate';
end

if ~strcmp(fromCat,toCat)
    error('mr.convert: fromUnit ''%s'' (%s) and toUnit ''%s'' (%s) are in different unit categories.', ...
        opt.fromUnit, fromCat, opt.toUnit, toCat);
end

% Convert to standard units
switch opt.fromUnit
    case {'Hz','Hz/m','Hz/m/s'}
        standard = in;
    case {'mT','mT/m'};
        standard = in*1e-3*opt.gamma;
    case {'uT'};
        standard = in*1e-6*opt.gamma;
    case 'rad/ms/mm'
        standard = in*1e6/(2*pi);
    case {'T','mT/m/ms','T/m/s'}
        standard = in*opt.gamma;
    case 'rad/ms/mm/ms'
        standard = in*1e9/(2*pi);
end

% Convert from standard units
switch opt.toUnit
    case {'Hz','Hz/m','Hz/m/s'}
        out = standard;
    case {'mT','mT/m'}
        out = 1e3*standard/opt.gamma;
    case 'uT'
        out = 1e6*standard/opt.gamma;
    case 'rad/ms/mm'
        out = standard*2*pi*1e-6;
    case {'T','mT/m/ms','T/m/s'}
        out = standard/opt.gamma;
    case 'rad/ms/mm/ms'
        out = standard*2*pi*1e-9;
end

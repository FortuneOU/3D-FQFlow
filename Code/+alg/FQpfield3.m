function [RP,param,SPECT,IDX] = FQpfield3(varargin)
%FQPFIELD3   3-D RMS acoustic pressure field of a planar 2-D transducer array
%   RP = FQpfield3(X,Y,Z,DELAYS,PARAM) computes the 3D radiation
%   pattern of a planar 2-D array excited with prescribed transmit delays.
%   The output is the **root-mean-square (RMS) acoustic pressure field**.
%
%   PARAM is a transducer/medium parameter structure.
%   The simulation is performed for observation points at coordinates (X,Y,Z).
%
%   Use PFIELD (2D version) for linear/convex arrays.
%
%   Units: X,Y,Z in meters; DELAYS in seconds.
%
%   DELAYS can be a matrix (M×N). This allows simulating multiple-plane
%   transmit wavefronts (MPT). Each row corresponds to one transmission law.
%
%   Notes:
%     - Use TXDELAY3 to generate delays for point/line/plane/diverging focus.
%     - FQpfield3 is called internally by SIMUS3 when simulating RF signals.
%
% -----------------------------------------------------------------------
% Coordinate Convention:
%   - X-axis is parallel to transducer surface, increasing left→right
%   - Y-axis is parallel to aperture plane, forming right-handed system
%   - Z-axis is perpendicular, positive in depth direction
% -----------------------------------------------------------------------
% Simplifications:
%   - Transducer elements directivity assumed frequency-independent by default
%   - Medium assumed homogeneous with scalar sound speed
%   - Valid for large distances where far-field approximations hold
% -----------------------------------------------------------------------


assert(nargin>0,'Do not have inputs1');
MyFolder = pwd;
narginchk(5,6)

%-- Input variables: X,Y,Z,DELAYS,PARAM,OPTIONS
x = varargin{1};
y = varargin{2};
z = varargin{3};
assert(isequal(size(x),size(y),size(z)),...
    'X, Y, and Z must be of same size.')
delaysTX = varargin{4};
param = varargin{5};
if nargin==6
    options = varargin{6};
else
    options = [];
end

%-- Check the transmit delays
assert(isnumeric(delaysTX) && all(delaysTX(~isnan(delaysTX))>=0),...
    'DELAYS must be a nonnegative array.')
if isvector(delaysTX) % we need a row vector
    delaysTX = delaysTX(:).'; 
end
NumberOfElements = size(delaysTX,2);

%--
% Note: delaysTX can be a matrix. This option can be used for MPT
% (multi-plane transmit) for example. In this case, each row represents a
% delay series. For example, for a 4-MPT sequence with a 1024-element
% matrix array, delaysTX has 4 rows and 1024 columns, i.e. size(delaysTX) =
% [4 1024].
%--
delaysTX = delaysTX.';

%-- Check if PFIELD3 is called by SIMUS3
isSIMUS3 = false;
if isfield(options,'CallFun')
    isSIMUS3 = strcmpi(options.CallFun,'simus3');
end

%-- Check if SIMUS3 is running on a parallel pool
try
    onppool = isSIMUS3 && options.ParPool;
catch
    onppool = false;
end



%---------------------------%
% Check the PARAM structure %
%---------------------------%

param = IgnoreCaseInFieldNames(param);

%-- 1) Center frequency (in Hz)
assert(isfield(param,'fc'),...
    'A center frequency value (PARAM.fc) is required.')
fc = param.fc; % central frequency (Hz)

%-- 2) Coordinates of the transducer elements (xe,ye)
% note: ze = 0 in this PFIELD3 version.
assert(isfield(param,'elements'),...
    ['PARAM.elements must contain the x- and y-locations ',...
    'of the transducer elements.'])
assert(size(param.elements,1)==2,...
    ['PARAM.elements must have two rows that contain the ',...
    'x (1st row) and y (2nd row) coordinates of the transducer elements.'])
xe = param.elements(1,:);
ye = param.elements(2,:);
assert(numel(xe)==NumberOfElements,...
    'The number of elements must match the number of transmit delays.')

%-- 3) Element width (in m)
assert(isfield(param,'width'),...
    'An element width (PARAM.width) is required.')
ElementWidth = param.width;
assert(isnumeric(ElementWidth) && isscalar(ElementWidth) &&...
    ElementWidth>0,'The element width must be positive.')

%-- 4) Element height (in m)
assert(isfield(param,'height'),...
    'An element height (PARAM.height) is required with PFIELD3 and SIMUS3.')
ElementHeight = param.height;
assert(isnumeric(ElementHeight) && isscalar(ElementHeight) &&...
    ElementHeight>0,'The element height must be positive.')

%-- 5) Fractional bandwidth at -6dB (in %)
if ~isfield(param,'bandwidth')
    param.bandwidth = 75;
end
assert(param.bandwidth>0 && param.bandwidth<200,...
    'The fractional bandwidth at -6 dB (PARAM.bandwidth, in %) must be in ]0,200[')

%-- 6) Baffle
%   An obliquity factor will be used if the baffle is not rigid
%   (default = SOFT baffle)
if ~isfield(param,'baffle')
    param.baffle = 'soft'; % default
end
if strcmpi(param.baffle,'rigid')
    NonRigidBaffle = false;
elseif strcmpi(param.baffle,'soft')
    NonRigidBaffle = true;
elseif isscalar(param.baffle)
    assert(param.baffle>0,...
        'The ''baffle'' field scalar must be positive')
    NonRigidBaffle = true;
else
    error('The ''baffle'' field must be ''rigid'',''soft'' or a positive scalar')
end

%-- 7) Longitudinal velocity (in m/s)
if ~isfield(param,'c')
    param.c = 1540; % default value
end
c = param.c; % speed of sound (m/s)

%-- 8) Attenuation coefficient (in dB/cm/MHz)
if ~isfield(param,'attenuation') % no attenuation, alpha_dB = 0
    param.attenuation = 0;
    alpha_dB = 0;
else
    alpha_dB = param.attenuation;
    assert(isscalar(alpha_dB) && isnumeric(alpha_dB) && alpha_dB>=0,...
        'PARAM.attenuation must be a nonnegative scalar')
end

%-- 9) Transmit apodization (no unit)
if ~isfield(param,'TXapodization')
    param.TXapodization = ones(1,NumberOfElements);
else
    assert(isvector(param.TXapodization) && isnumeric(param.TXapodization),...
        'PARAM.TXapodization must be a vector')
    assert(numel(param.TXapodization)==NumberOfElements,...
        'PARAM.TXapodization must be of length = (number of elements)')
end
% apodization is 0 where TX delays are NaN:
idx = isnan(delaysTX);
param.TXapodization(any(idx,2)) = 0;
delaysTX(idx) = 0;

%-- 10) TX pulse: Number of wavelengths
if ~isfield(param,'TXnow')
    param.TXnow = 1;
end
NoW = param.TXnow;
assert(isscalar(NoW) && isnumeric(NoW) && NoW>0,...
    'PARAM.TXnow must be a positive scalar.')

%-- 11) TX pulse: Frequency sweep for a linear chirp
if ~isfield(param,'TXfreqsweep') || isinf(NoW)
    param.TXfreqsweep = [];
end
FreqSweep = param.TXfreqsweep;
assert(isempty(FreqSweep) ||...
    (isscalar(FreqSweep) && isnumeric(FreqSweep) && FreqSweep>0),...
    'PARAM.TXfreqsweep must be empty (windowed sine) or a positive scalar (linear chirp).')

%----------------------------------%
% END of Check the PARAM structure %
%----------------------------------%



%-----------------------------%
% Check the OPTIONS structure %
%-----------------------------%

options = IgnoreCaseInFieldNames(options);

%-- 1) dB threshold
%     (in dB: faster computation if lower value, but less accurate)
if ~isfield(options,'dBThresh')
    options.dBThresh = -60; % default is -60dB in PFIELD3
end
assert(isscalar(options.dBThresh) && isnumeric(options.dBThresh) &&...
    options.dBThresh<=0,'OPTIONS.dBThresh must be a nonpositive scalar.')

%-- 2) Frequency-dependent directivity?
if isfield(options,'FullFrequencyDirectivity')
    isFFD = options.FullFrequencyDirectivity;
else
    isFFD = false; % default
    % By default, the directivity of the elements depends on the center
    % frequency only. This makes the algorithm faster. 
end
assert(isscalar(isFFD) && islogical(isFFD),...
    'OPTIONS.FullFrequencyDirectivity must be a logical scalar (true or false).')

%-- 3) Element splitting
%
% --- A short note about the algorithm:
% Far-field equations are used in PFIELD3. Each transducer element of the
% array is split into M-by-N small rectangles, so that these MN rectangles
% have size smaller than one wavelength by one wavelength. The far-field
% condition is acceptable for these small rectangles.
%---
if isfield(options,'ElementSplitting') && ~isempty(options.ElementSplitting)
    assert(numel(options.ElementSplitting)==2,...
        'OPTIONS.ElementSplitting must be a two-element vector.')
    M = options.ElementSplitting(1);
    N = options.ElementSplitting(2);
    assert(isscalar(M) & M==round(M) & M>0 &...
        isscalar(N) & N==round(N) & N>0,...
        'OPTIONS.ElementSplitting must contain two positive integers.')
else
    LambdaMin = c/(fc*(1+param.bandwidth/200));
    M = ceil(ElementWidth/LambdaMin);
    N = ceil(ElementHeight/LambdaMin);
end

%-- 4) Wait bar
if ~isfield(options,'WaitBar')
    options.WaitBar = true;
end
assert(isscalar(options.WaitBar) && islogical(options.WaitBar),...
    'OPTIONS.WaitBar must be a logical scalar (true or false).')

%-- Advanced (masked) options: Frequency step (scaling factor)
% The frequency step is determined automatically. It is tuned to avoid
% significant interferences due to unadapted discretization. The frequency
% step can also be adjusted by using a scaling factor. For a fast check,
% you may use a scaling factor>1. For a smoother result, you may use a
% scaling factor<1.
if ~isfield(options,'FrequencyStep')
    options.FrequencyStep = 1;
end
assert(isscalar(options.FrequencyStep) &&...
    isnumeric(options.FrequencyStep) && options.FrequencyStep>0,...
    'OPTIONS.FrequencyStep must be a positive scalar.')

%------------------------------------%
% END of Check the OPTIONS structure %
%------------------------------------%


%-----
% SIMUS3 first runs PFIELD3 with empty X,Y,Z to detect syntax errors.
if isSIMUS3 && isempty(x), RP = []; return, end
%-----


%------------------------------------%
% POINT LOCATIONS, DISTANCES & GRIDS %
%------------------------------------%

siz0 = size(x);
nx = numel(x);

%-- Coordinates of the points where pressure is needed
x = x(:); y = y(:); z = z(:);

%-- Cast x, y, and z to single class
x = cast(x,'single');
y = cast(y,'single');
z = cast(z,'single');

%-- Centroids of the sub-elements
%-- note: Each elements is split into M-by-N sub-elements.
% X-position (xi) and Y-position (yi) of the centroids of the sub-elements
% (relative to the centers of the transducer elements).
% The values in xi are in the range ]-ElementWidth/2 ElementWidth/2[.
% The values in yi are in the range ]-ElementHeight/2 ElementHeight/2[.
% (if M = 1 and N = 1, then xi = yi = 0).
SegWidth = ElementWidth/M;
xi = -ElementWidth/2 + SegWidth/2 + (0:M-1)*SegWidth;
SegHeight = ElementHeight/N;
yi = -ElementHeight/2 + SegHeight/2 + (0:N-1)*SegHeight;
[xi,yi] = meshgrid(xi,yi);
xi = reshape(xi,[1 1 M*N]);
yi = reshape(yi,[1 1 M*N]);

%-- Out-of-field points
% Null pressure will be assigned to out-of-field points.
isOUT = z<0;

%-- Variables that we need:
%
% Note: We work in an ISO spherical system for each sub-element
%       r = distance between the segment centroid and the point of interest
%       sinT = sine theta: theta is the polar angle.
%       cosT = cosine theta.
%       sinP = sine phi: phi is the azimuthal angle.
%       cosP = cosine phi.
%       They are of size [numel(x) NumberOfElements M*N].
%
dxi = x-xi-xe;
dyi = y-yi-ye;
d2 = dxi.^2+dyi.^2;
r = sqrt(d2+z.^2);
%---
epss = eps('single');
cosT = (z+epss)./(r+epss);
sinT = (sqrt(d2)+epss)./(r+epss);
cosP = (dxi+epss)./(sqrt(d2)+epss);
sinP = (dyi+epss)./(sqrt(d2)+epss);
clear dxi dyi d2
%---
% The term 1/r is present in the equations (problems if r is very small!):
% small r values are replaced by lambda/2
lambda = c/fc;
r(r<lambda/2) = lambda/2;

%-------------------------------------------%
% end of POINT LOCATIONS, DISTANCES & GRIDS %
%-------------------------------------------%


mysinc = @(x) sin(abs(x)+epss)./(abs(x)+epss); % cardinal sine
% [note: In MATLAB, sinc is sin(pi*x)/(pi*x)]

%-------------------%
% FREQUENCY SPECTRA %
%-------------------%

%-- FREQUENCY SPECTRUM of the transmitted pulse
if isempty(FreqSweep)
    % We want a windowed sine of width PARAM.TXnow
    T = NoW/fc; % temporal pulse width
    if isinf(T); T = 1e6; end
    wc = 2*pi*fc;
    pulseSpectrum = @(w) 1i*(mysinc(T*(w-wc)/2)-mysinc(T*(w+wc)/2));
else
    % We want a linear chirp of width PARAM.TXnow
    % (https://en.wikipedia.org/wiki/Chirp_spectrum#Linear_chirp)
    T = NoW/fc; % temporal pulse width
    if isinf(T); T = 1e6; end
    wc = 2*pi*fc;
    dw = 2*pi*FreqSweep;
    s2 = @(w) sqrt(pi*T/dw)*exp(-1i*(w-wc).^2*T/2/dw).*...
        (fresnelint((dw/2+w-wc)/sqrt(pi*dw/T)) + fresnelint((dw/2-w+wc)/sqrt(pi*dw/T)));
    pulseSpectrum = @(w) (1i*s2(w)-1i*s2(-w))/T;
end

%-- FREQUENCY RESPONSE of the ensemble PZT + probe
% We want a generalized normal window (6dB-bandwidth = PARAM.bandwidth)
% (https://en.wikipedia.org/wiki/Window_function#Generalized_normal_window)
wB = param.bandwidth*wc/100; % angular frequency bandwidth
p = log(126)/log(2*wc/wB); % p adjusts the shape
probeSpectrum = @(w) exp(-(abs(w-wc)/(wB/2/log(2)^(1/p))).^p);
% The frequency response is a pulse-echo (transmit + receive) response. A
% square root is thus required when calculating the pressure field:
probeSpectrum = @(w) sqrt(probeSpectrum(w));
% Note: The spectrum of the pulse (pulseSpectrum) will be then multiplied
% by the frequency-domain tapering window of the transducer (probeSpectrum)

%-- FREQUENCY STEP
if isSIMUS3 % PFIELD3 has been called by SIMUS3
    df = options.FrequencyStep;
else % We are in PFIELD3 only (i.e. not called by SIMUS3)
    % The frequency step df is chosen to avoid interferences due to
    % inadequate discretization.
    % -- df = frequency step (must be sufficiently small):
    % One has exp[-i(k r + w delay)] = exp[-2i pi(f r/c + f delay)] in the Eq.
    % One wants: the phase increment 2pi(df r/c + df delay) be < 2pi.
    % Therefore: df < 1/(r/c + delay).
    df = 1/(max(r(:)/c) + max(delaysTX(:)));
    df = options.FrequencyStep*df;
    % note: df is here an upper bound; it will be recalculated below
end

%-- FREQUENCY SAMPLES
Nf = 2*ceil(param.fc/df)+1; % number of frequency samples
f = linspace(0,2*param.fc,Nf); % frequency samples
df = f(2); % update the frequency step
%- we keep the significant components only by using options.dBThresh
S = abs(pulseSpectrum(2*pi*f).*probeSpectrum(2*pi*f));
GdB = 20*log10(S/max(S)); % gain in dB
IDX = GdB>options.dBThresh;
IDX(find(IDX,1):find(IDX,1,'last')) = true;
f = f(IDX);
nSampling = length(f);

%-- we need VECTORS
pulseSPECT = pulseSpectrum(2*pi*f); % pulse spectrum
probeSPECT = probeSpectrum(2*pi*f); % probe response

%--------------------------%
% end of FREQUENCY SPECTRA %
%--------------------------%



%-- Wait bar
options.WaitBar = options.WaitBar & (nSampling>10);
if options.WaitBar
    if isSIMUS3
        wbtitle = 'Let SIMUS3 do the work for you...';
        wbname = 'SIMUS3 / www.biomecardio.com';
    else
        wbtitle = 'Let PFIELD3 do the work for you...';
        wbname = 'PFIELD3 / www.biomecardio.com';
    end
    hwb = waitbar(0,wbtitle,'Name',wbname);
end

%-- Initialization
RP = 0; % RP = Radiation Pattern
if isSIMUS3
    %- For SIMUS3 only (we need the full spectrum of RX signals):
    SPECT = zeros([nSampling NumberOfElements],'like',single(1i));
elseif nargout==3
    SPECT = zeros([nSampling nx],'like',single(1i));
end

%-- Obliquity factor (baffle property)
%   An obliquity factor is required if the baffle is not rigid.
%   [Th = angle relative to the element normal axis]
if NonRigidBaffle
    if strcmpi(param.baffle,'soft')
        ObliFac = cosT;
    else % param.baffle is a scalar
        ObliFac = cosT./(cosT+param.baffle);
    end
else % 1 if rigid baffle
    ObliFac = ones(size(cosT));
end

%-- Note on Attenuation
% Reference: Diagnostic ultrasound imaging - inside out (T.L. Szabo)
%            Chapter 4: Attenuation
% Key reference: Acoustics for ultrasound imaging (Ben Cox, 2013)
%                Chapter 5: Acoustic attenuation and absorption
% We will use this attenuation-based wavenumber:
%   kwa = alpha_dB/8.69*f(k)/1e6*1e2; % P(z,f) = P0 exp(-alpha*f*z/8.69)
%   note: 20/log(10) ~ 8.69

%-- EXPONENTIAL arrays of size [numel(x) NumberOfElements M]
kw = 2*pi*f(1)/c; % wavenumber
kwa = alpha_dB/8.69*f(1)/1e6*1e2; % attenuation-based wavenumber
EXP = exp(-kwa*r + 1i*mod(kw*r,2*pi)); % faster than exp(-kwa*r+1i*kw*r)
%-- Exponential array for the increment wavenumber dk
dkw = 2*pi*df/c;
dkwa = alpha_dB/8.69*df/1e6*1e2;
EXPdf = exp((-dkwa + 1i*dkw)*r);
%lnEXPdf = (-dkwa + 1i*dkw)*r;

%-- We replace EXP by EXP.*ObliFac./r
EXP = EXP.*ObliFac./r;
%lnEXP = -kwa*r + 1i*mod(kw*r,2*pi) + log(ObliFac./r); 
clear ObliFac r
gpu = gpuDevice;
disp(['GPU All  ',num2str(gpu.AvailableMemory/1e6),'MB memory']);
EXP_gpu = gpuArray(single(EXP));
EXPdf_gpu = gpuArray(single(EXPdf));




disp(['GPU remain  ',num2str(gpu.AvailableMemory/1e6),'MB memory']);


%-- TX apodization
APOD = param.TXapodization(:);

%-- Simplified directivity (if not dependent on frequency)
% In the "simplified directivity" version, the directivities of the
% sub-elements depend on the center frequency ONLY. It is thus not needed
% to calculate the directivity arrays (DIRx and DIRy) in the following
% for-loop. These directivities DIRx and DIRy are included in the variable
% EXP to reduce storage.
if ~isFFD
    kc = 2*pi*fc/c; % center wavenumber
    DIRx = mysinc(kc*SegWidth/2*cosP.*sinT); % x-directivity of each segment
    DIRy = mysinc(kc*SegHeight/2*sinP.*sinT); % y-directivity of each segment
    EXP = EXP.*DIRx.*DIRy;
    clear DIRx DIRy
end
 


%-----------------------------%
% SUMMATION OVER THE SPECTRUM %
%-----------------------------%

tstart = tic;
for k = 1:nSampling

    kw = 2*pi*f(k)/c; % wavenumber
    
    %-- Exponential array of size [numel(x) NumberOfElements MxN]
    % For all k, we need: EXP = exp((-kwa+1i*kw)*r)
    %                         with kw = 2*pi*f(k)/c;
    %                     and with kwa = alpha_dB/8.7*f(k)/1e6*1e2;
    % Since f(k) = f(1)+(k-1)df, we use the following recursive product:
    if k>1
        %EXP = EXP.*EXPdf;

      %  disp(['GPU1_1 ---remain  ',num2str(gpu.AvailableMemory/1e6),'MB memory']);
        EXP_gpu = EXP_gpu .* EXPdf_gpu;
       % disp(['GPU1_2 ---remain  ',num2str(gpu.AvailableMemory/1e6),'MB memory']);

        %disp('here');
    end
        
    %-- Directivity (if frequency-dependent)
    if isFFD % isFFD = true -> frequency-dependent directivity
        DIRx = mysinc(kw*SegWidth/2*cosP.*sinT); % x-directivity
        DIRy = mysinc(kw*SegHeight/2*sinP.*sinT); % y-directivity
        DIR = DIRx.*DIRy;
    end
        
    %-- Radiation patterns of the single elements
    % They are the combination of the far-field patterns of the M small
    % segments that make up the single elements
    %--
    if isFFD % isFFD = true -> frequency-dependent directivity
        RPmono = mean(DIR.*EXP,3); % summation over the M*N small segments
    else % isFFD = false: the directivity depends on center frequency only
         % note: the directivity (DIR) has already been included in EXP
        if M*N>1           
          
           %   RPmono = mean(EXP,3); % summation over the M*N small segments
            RPmono = mean(EXP_gpu,3); % summation over the M*N small segments
           

        else
            RPmono = EXP;
        end
    end
        
    %-- Transmit delays + Transmit apodization
    % use of SUM: summation over the number of delay series (e.g. MLT)
    DELAPOD = sum(exp(1i*kw*c*delaysTX),2).*APOD;
    
    %-- Summing the radiation patterns generating by all the elements
    RPk = RPmono*DELAPOD;
    %- include spectrum responses:
    RPk = pulseSPECT(k)*probeSPECT(k)*RPk;
    RPk(isOUT) = 0;
       
    %-- Output
    if isSIMUS3 % Receive: for SIMUS3 only (spectra of the RF signals)
        SPECT(k,:) = probeSPECT(k) *... % the array bandwidth is considered
            (RPk.*options.RC(:)).'*RPmono ... % pressure received by the elements
            ; % *f(k)^2/fc^2; % Rayleigh scattering (OPTIONAL)
        if any(param.RXdelay) % reception delays, if any
            SPECT(k,:) = SPECT(k,:).*exp(1i*kw*c*param.RXdelay);
        end
    else % using PFIELD3 alone
        RP = RP + abs(RPk).^2; % acoustic intensity
        if nargout==3
            SPECT(k,:) = RPk;
        end
    end
    
    %- update the wait bar
    if options.WaitBar && rem(k,10)==0
        tstep = toc(tstart);
        trem = tstep*(nSampling/k-1);
        waitbar(k/nSampling,hwb,...
            {['Elapsed: ',int2str(floor(tstep/60)) ,' min ',...
            int2str(floor(rem(tstep,60))),' s'],...
            ['Remaining: ',int2str(floor(trem/60)) ,' min ',...
            int2str(floor(rem(trem,60))),' s']})
    end
    
end

%------------------------------------%
% end of SUMMATION OVER THE SPECTRUM %
%------------------------------------%



% Close the wait bar
if options.WaitBar, close(hwb), end

% Correcting factor (including integration step, df)
if isinf(NoW)
    CorFac = 1;
else
    CorFac = df;
end
% if exist('SPECT','var'), SPECT = SPECT*CorFac; end
RP = RP*CorFac;

% RMS acoustic pressure (if we are in PFIELD3 only)
if ~isSIMUS3
    RP = reshape(sqrt(RP),siz0);
    if nargout>2
        SPECT = reshape(SPECT.',[siz0 nSampling]);
    end
end

end







function f = fresnelint(x)

% FRESNELINT Fresnel integral.
%
% J = FRESNELINT(X) returns the Fresnel integral J = C + 1i*S.
%
% We use the approximation introduced by Mielenz in
%       Klaus D. Mielenz, Computation of Fresnel Integrals. II
%       J. Res. Natl. Inst. Stand. Technol. 105, 589 (2000), pp 589-590
%

siz0 = size(x);
x = x(:);

issmall = abs(x)<=1.6;
c = zeros(size(x));
s = zeros(size(x));

% When |x| < 1.6, a Taylor series is used (see Mielenz's paper)
if any(issmall)
    n = 0:10;
    cn = [1 cumprod(-pi^2*(4*n+1)./(4*(2*n+1).*(2*n+2).*(4*n+5)))];
    sn = [1 cumprod(-pi^2*(4*n+3)./(4*(2*n+2).*(2*n+3).*(4*n+7)))]*pi/6;
    n = [n 11];
    c(issmall) = sum(cn.*x(issmall).^(4*n+1),2);
    s(issmall) = sum(sn.*x(issmall).^(4*n+3),2);    
end

% When |x| > 1.6, we use the following:
if any(~issmall)
    n = 0:11;
    fn = [0.318309844, 9.34626e-8, -0.09676631, 0.000606222, ...
        0.325539361, 0.325206461, -7.450551455, 32.20380908, ...
        -78.8035274, 118.5343352, -102.4339798, 39.06207702];
    gn = [0, 0.101321519, -4.07292e-5, -0.152068115, -0.046292605, ...
        1.622793598, -5.199186089, 7.477942354, -0.695291507, ...
        -15.10996796, 22.28401942, -10.89968491];
    fx = sum(fn.*x(~issmall).^(-2*n-1),2);
    gx = sum(gn.*x(~issmall).^(-2*n-1),2);    
    c(~issmall) = 0.5*sign(x(~issmall)) + ...
        fx.*sin(pi/2*x(~issmall).^2) - gx.*cos(pi/2*x(~issmall).^2);
    s(~issmall) = 0.5*sign(x(~issmall)) - ...
        fx.*cos(pi/2*x(~issmall).^2) - gx.*sin(pi/2*x(~issmall).^2);
end

f = reshape(c,siz0) + 1i*reshape(s,siz0);

end


function structArray = IgnoreCaseInFieldNames(structArray)

switch inputname(1)
    case 'param'
        fieldLIST = {'attenuation','baffle','bandwidth','c','elements','fc',...
            'fnumber','focus','fs','height','kerf','movie','Nelements',...
            'passive','pitch','radius','RXangle','RXdelay',...
            'TXapodization','TXfreqsweep','TXnow','t0','width'};
    case 'options'
        if isstruct(structArray)
            fieldLIST = {'dBThresh','ElementSplitting',...
                'FullFrequencyDirectivity','FrequencyStep','ParPool',...
                'WaitBar'};
        else
            return
        end
end

OldFieldNames = fieldnames(structArray);
tmp = lower(OldFieldNames);
assert(length(tmp)==length(unique(tmp)),...
    ['The structure ' upper(inputname(1)),...
    ' contains duplicate field names (when ignoring case).'])

[idx,loc] = ismember(lower(fieldLIST),tmp);
idx = find(idx); loc = loc(idx);
for k = 1:length(idx)
    tmp = eval(['structArray.' OldFieldNames{loc(k)}]); %#ok
    structArray = rmfield(structArray,OldFieldNames{loc(k)});
    eval(['structArray.' fieldLIST{idx(k)} ' = tmp;']) %#ok
end

end


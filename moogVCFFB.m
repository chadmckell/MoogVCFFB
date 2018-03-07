%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Moog VCF Filter (Forward & Backward Approximations)
% 
% Author: Chad McKell
% Date: 27.03.17
%
% Description: Virtual analog model of the Moog VCF ladder filter. This
% script uses forward and backward integrators to approximate the impulse
% response of the filter. The approximated transfer functions are then
% compared with the exact calculation of the Moog VCF transfer function.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
tic; close; clc; clear;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define global variables
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

SR = 44100; % sample rate [samples/sec]
k = 1/SR; % time step [sec]
f0 = 1000; % resonant frequency [Hz]
w0 = 2*pi*f0; % angular frequency [rad/sec]
r = 0.5; % tuning parameter (a number between 0 and 1) 
Tf = 3; % total time length of simulation [sec]
Nf = floor(Tf*SR); % total sample length of simulation [samples]
A = w0*[-1 0 0 -4*r; 1 -1 0 0; 0 1 -1 0; 0 0 1 -1]; % system matrix
b = w0*[1; 0; 0; 0]; % forcing vector
c = [0; 0; 0; 1]; % state mixture vector

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Error handling: Terminate code for the following errors
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% r is out of bounds
if r > 1 || r < 0
   error('r must be a number in the range [0,1].')
end

% SR has a non-zero decimal
if rem(SR,1) ~= 0
   error('SR must be zero-decimal.')
end

% SR, f0, Tf, or r is not a number
if ~isfloat(SR) || ~isfloat(f0) || ~isfloat(Tf) || ~isfloat(r) 
   error('SR, f0, Tf, and r must each be a number.')
end

% SR, f0, Tf, or r is not real
if ~isreal(SR) || ~isreal(f0) || ~isreal(Tf) || ~isreal(r)   
   error('SR, f0, Tf, and r must each be real.')
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%  Initialize the input sequence
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
u = zeros(Nf,1); 
u(1) = 1;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%  Implement forward integration method
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Initialize state vectors
xForward = zeros(4,1); % current time step
xForward1 = zeros(4,1); % one time step back

% Initialize output vector
yForward = zeros(Nf,1);

% Set coefficient
coef = eye(4) + k*A;

% Stability check: k is out of bounds
for n=0:3 
    
    % Compute real value of each eigen value of A
    eigen = real(w0*(-1 + sqrt(2)*r^(1/4)*exp(1i*(pi/4 + n*pi/2))));
    if eigen > 0
       error('The real part of each eigen value of A must be negative.')
    
    % Calculate stability condition
    else
        if k > -2/eigen
           error(['k cannot be greater than ' num2str(-2/eigen) '.'])
        end
    end
end

for n = 2:Nf
   
   % Main algorithm 
   xForward = coef*xForward1 + k*b*u(n-1);
   yForward(n) = c'*xForward;
   
   % Set value of xForward1 equal to next grid line
   xForward1 = xForward;
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%  Implement backward integration method
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Initialize state vectors
xBackward = zeros(4,1); % current time step
xBackward1 = zeros(4,1); % one time step back

% Initialize output vector
yBackward = zeros(Nf,1);

% Compute coefficient Q
Q = eye(4) - k*A;

for n = 1:Nf
   
   % Main algorithm
   d = xBackward1 + k*b*u(n);
   xBackward = Q\d;
   yBackward(n) = c'*xBackward;
   
   % Set value of xBackward1 equal to next grid line
   xBackward1 = xBackward;
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%  Compute transfer functions of forward and backward output signals
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Calculate Fourier transforms
U = fft(u);
YForward = fft(yForward);
YBackward = fft(yBackward);

% Compute transfer functions
HForward = U./YForward;
HBackward = U./YBackward;

% Convert to loglog form
HForwardLog = log(HForward) - max(log(HForward));
HBackwardLog = log(HBackward) - max(log(HBackward));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%  Compute exact transfer function
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

fmax = 22500;
fbin = 0:fmax/1000:fmax;
fnum = length(fbin);
Hexact = zeros(fnum,1);

for n = 1:fnum
    inv = (2*pi*1i*fbin(n)*eye(4)-A)\b;
    Hexact(n) = c'*inv;
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Define plotting parameters 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

n = 0:Nf-1; % length bins [samples]
t = (n/SR)'; % time bins [sec]
fk = n'*SR/Nf; % frequency bins [Hz]
nyquist = SR/2; % Nyquist frequency [Hz]
font = 14; % font size for plots

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Generate plot of transfer functions
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

figure('units','normalized','outerposition',[0 0 1 1]);
orient landscape
plot(fk, abs(HBackwardLog), 'LineWidth', 2);
hold on
plot(fk, abs(HForwardLog), 'LineWidth', 2);
hold on
plot(fbin, real(log(Hexact)), 'LineWidth', 2);
xlim([0 nyquist]);
title('Transfer Functions of Moog VCF Ladder Filter');
xlabel('Frequency (Hz)'); 
ylabel('Magnitude');
legend('Forward', 'Backward', 'Exact')
set(gca,'fontsize',font)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Check code efficiency
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
toc % print elapsed time



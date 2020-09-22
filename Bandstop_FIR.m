f_samp = 1200e3;

%Band Edge speifications
fp1 = 258e3;
fs1 = 278e3;
fs2 = 333e3;
fp2 = 353e3;
t_w = 10e3;
del_w_t = t_w*2*pi/f_samp;
Wc1 = fp1*2*pi/f_samp;
Wc2  = fp2*2*pi/f_samp;
%Kaiser paramters
A = -20*log10(0.15);
beta = 0;

%Wn = [(fs1+fp1)/2 (fs2+fp2)/2]*2/f_samp;        %average value of the two paramters
N_min = ceil((A-7.95) / (2*2.285*del_w_t));       %empirical formula for N_min

%Window length for Kaiser Window
n=N_min + 13

%% Ideal bandstop impulse response of length "n"

bs_ideal =  ideal_lp(pi,n) -ideal_lp(Wc2 - del_w_t,n) + ideal_lp(Wc1 + del_w_t,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandStop = bs_ideal .* kaiser_win;
fvtool(FIR_BandStop);         %frequency response

%% magnitude response
[H,f] = freqz(FIR_BandStop,1,f_samp/2, f_samp);
figure(1), plot(f,abs(H), 'LineWidth', 3)
hold on
plot(fp1,abs(H(fp1)),'r*')
plot(fp2,abs(H(fp2)),'r*')
plot(fs1,abs(H(fs1)),'g*')
plot(fs2,abs(H(fs2)),'g*')
grid
ylabel('Magnitude');
xlabel('Frequency (Hz)');
title('Magnitude Plot');
disp(n)
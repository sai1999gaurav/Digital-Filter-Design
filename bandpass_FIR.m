f_samp = 1200e3;

%Band Edge speifications
fs1 = 251e3;
fp1 = 271e3;
fp2 = 316e3;
fs2 = 336e3;
t_w = 10e3;
del_w_t = t_w*2*pi/f_samp;
Wc1 = fp1*2*pi/f_samp;
Wc2  = fp2*2*pi/f_samp;
%Kaiser paramters
A = -20*log10(0.15);
beta = 0;
% if(A < 21)
%     beta = 0;
% elseif(A <51)
%     beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
% else
%     beta = 0.1102*(A-8.7);
% end

N_min = ceil((A-7.95) / (2*2.285*del_w_t));           %empirical formula for N_min
n = N_min + 17;
% Window length for Kaiser Window
%% 
%while(1)

%Ideal bandpass impulse response of length "n"
bp_ideal = ideal_lp(Wc2 + del_w_t,n) - ideal_lp(Wc1 - del_w_t,n);
%disp(bp_ideal)
%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandPass = bp_ideal .* kaiser_win;
fvtool(FIR_BandPass);         %frequency response

%% magnitude response
[H,f] = freqz(FIR_BandPass,1,f_samp/2, f_samp);
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

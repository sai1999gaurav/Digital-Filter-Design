%Butterworth Analog LPF parameters
Wc = 1.1;              %cut-off frequency
N = 7;                  %order 

%poles of Butterworth polynomial of degree 7 in the open CLHP from wolfram 
p1 = -1.1;
p2 = -0.991066 - 0.477272i;
p3 = -0.991066 + 0.477272i;
p4 = -0.685839 - 0.860015i;
p5 = -0.685839 + 0.860015i;
p6 = -0.244773 - 1.072420i;
p7 = -0.244773 + 1.072420i;
%Band Edge speifications in kHz
fs1 = 251;
fp1 = 271;
fp2 = 316;
fs2 = 336;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 1200;         
ws1 = tan(fs1*pi/f_samp);
wp1 = tan(fp1*pi/f_samp);
wp2 = tan(fp2*pi/f_samp);
ws2 = tan(fs2*pi/f_samp);

%Parameters for Bandpass Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

[num,den] = zp2tf([],[p1 p2 p3 p4 p5 p6 p7],Wc^N);    %TF with poles p1-p8 and numerator Wc^N and no zeroes
                                                        %numerator chosen to make the DC Gain = 1
disp(num)
disp(den)
%% Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bpf(s) = analog_lpf((s*s + W0*W0)/(B*s));        %bandstop transformation
discrete_bpf(z) = analog_bpf((z-1)/(z+1));              %bilinear transformation

%% coeffs of analog bpf
[ns, ds] = numden(analog_bpf(s));                   %numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;
disp(ns);
disp(ds);

%% coeffs of discrete bpf
[nz, dz] = numden(discrete_bpf(z));                     %numerical simplification to collect coeffs                    
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %coeffs to matrix form
k = dz(1);                                              %normalisation factor
dz = dz/k;
nz = nz/k;
disp('Nz');
disp(nz);
disp(dz);
fvtool(nz,dz)                                           %frequency response

%% magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,600000, f_samp*1e3);
plot(f,abs(H), 'LineWidth', 3)
hold on
plot(fp1*1e3,abs(H(fp1*1e3)),'r*')
x = fp1*1e3;
y = abs(H(fp1*1e3));
textString = sprintf(' (%d, %.2f)', fp1, y);
text(x,y,textString, 'FontSize',9.5)
plot(fp2*1e3,abs(H(fp2*1e3)),'r*')
x = fp2*1e3;
y = abs(H(fp2*1e3));
textString = sprintf(' (%d, %.2f)', fp2, y);
text(x,y,textString, 'FontSize',9.5)
plot(fs1*1e3,abs(H(fs1*1e3)),'g*')
x = fs1*1e3;
y = abs(H(fs1*1e3));
textString = sprintf(' (%d, %.2f)', fs1, y);
text(x,y,textString, 'FontSize',9.5)
plot(fs2*1e3,abs(H(fs2*1e3)),'g*')
x = fs2*1e3;
y = abs(H(fs2*1e3));
textString = sprintf(' (%d, %.2f)', fs2, y);
text(x,y,textString, 'FontSize',9.5)
grid
ylabel('Magnitude');
xlabel('Frequency (Hz)');
title('Magnitude Plot');
%legend('Magnitude Response', 'Pass band cut-off','','Stop band cut-off')
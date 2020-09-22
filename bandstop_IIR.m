%Chebyshev LPF parameters
D1 = 1/(0.85*0.85)-1;       %since delta is 0.15
epsilon = sqrt(D1);         %epsilon was set to this value to satisfy required inequality
N = 3;


%poles of Butterworth polynomial of degree 8 in the open CLHP 
p1 = -0.43105;
p2 = -0.21553 + 0.94306i;
p3 = -0.21553 - 0.94306i;

%evaluating the Transfer function of Chebyshev Analog LPF


%Band Edge speifications
fp1 = 258e3;
fs1 = 278e3;
fs2 = 333e3;
fp2 = 353e3;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 1200e3;    
wp1 = tan(fp1*pi/f_samp);
ws1 = tan(fs1*pi/f_samp);
ws2 = tan(fs2*pi/f_samp);
wp2 = tan(fp2*pi/f_samp);

%Parameters for Bandstop Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

[num,den] = zp2tf([],[p1 p2 p3],-p1*p2*p3);   %TF with poles p1-p8 and numerator Wc^N and no zeroes
                                                        %numerator chosen to make the DC Gain = 1
%% Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bsf(s) = analog_lpf((B*s)/(s*s + W0*W0));        %bandstop transformation
discrete_bsf(z) = analog_bsf((z-1)/(z+1));              %bilinear transformation

%coeffs of analog bsf
[ns, ds] = numden(analog_bsf(s));                   %numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k
ns = ns/k

%% coeffs of discrete bsf
[nz, dz] = numden(discrete_bsf(z));                     %numerical simplification to collect coeffs                    
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %coeffs to matrix form
k = dz(1);                                              %normalisation factor
dz = dz/k
nz = nz/(k)
fvtool(nz,dz)                                           %frequency response

%% magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,f_samp, f_samp);
plot(f,abs(H), 'LineWidth', 3)
hold on
ylabel('Magnitude');
xlabel('Frequency (Hz)');
title('Magnitude Plot');
grid
%Elliptic LPF parameters
N = 3;

%Band Edge speifications
fp1 = 258;
fs1 = 278;
fs2 = 333;
fp2 = 353;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 1200;    
wp1 = tan(fp1*pi/f_samp);
ws1 = tan(fs1*pi/f_samp);
ws2 = tan(fs2*pi/f_samp);
wp2 = tan(fp2*pi/f_samp);

%Parameters for Bandstop Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

num = 0.09*[1 0 2.2352]
den = [1 0.6436 0.6503 0.2012]
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
nz = nz/(k);
disp(nz/nz(1))
fvtool(nz,dz)                                           %frequency response

%% magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,600000, f_samp*1e3);
fp1 = 258;
fs1 = 278;
fs2 = 333;
fp2 = 353;
plot(f,abs(H), 'LineWidth', 3)
hold on
plot((fp1-15)*1e3,abs(H((fp1-15)*1e3)),'r*')
x = fp1*1e3;
y = abs(H((fp1-15)*1e3));
textString = sprintf(' (%d, %.2f)', fp1, y);
text(x,y,textString, 'FontSize',9.5)
plot((fp2+15)*1e3,abs(H((fp2+15)*1e3)),'r*')
x = fp2*1e3;
y = abs(H((fp2+15)*1e3));
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
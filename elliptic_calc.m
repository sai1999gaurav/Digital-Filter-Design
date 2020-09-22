%% 
Ap = 1.41;
n = 3;
Lambda = log( (10^(0.05*Ap) + 1)/(10^(0.05*Ap) - 1) )/(2*n)
q = 0.0245;
%num_series
n_s = 0;
for m = 0:3
   n_s = n_s + (-1.^m)*(q.^(m*(m+1)))*(sinh((2*m + 1)*Lambda));
end
n_s = n_s*2*(q.^0.25);
d_s = 1;
for m = 1:4
   d_s = d_s + 2*(-1.^m)*(q.^(m*m))*(cosh((2*m)*Lambda));
end
sigma_0 = abs(n_s/d_s)
%%
k = 0.57;
n = 3;
W = sqrt((1 + k*(sigma_0.^2))*(1 + (sigma_0.^2)/k))
%%
r = (n-1)/2;
omg = zeros(1,r);
for i = 1:r
    n_s = 0;
    for m = 0:3
        n_s = n_s + (-1.^m)*(q.^(m*(m+1)))*(sin((2*m + 1)*pi*i/n));
    end
    n_s = n_s*2*(q.^0.25);
    d_s = 1;    
    for m = 1:4
        d_s = d_s + 2*(-1.^m)*(q.^(m*m))*(cos((2*m)*pi*i/n));
    end 
    omg(1,i) = n_s./d_s    
end
%%
V = zeros(1,r);
for i = 1:r
    V(1,i) = sqrt((1 - k*(omg(1,i).^2))*(1 - (omg(1,i).^2)/k))
end
%%
a0 = zeros(1,r);
b0 = zeros(1,r);
b1 = zeros(1,r);
for i = 1:r
    a0(1,i) = 1/(omg(1,i).^2)
    b0(1,i) = ( (sigma_0*V(1,i)).^2 + (omg(1,i)*W).^2 )/((1 + (sigma_0.^2)*(omg(1,i).^2)).^2)
    b1(1,i) = (2*sigma_0*V(1,i))/(1 + (sigma_0*omg(1,i)).^2)
end
%%
H0 = sigma_0;
for i = 1:r
    H0 = H0*b0(1,i)/a0(1,i)    
end
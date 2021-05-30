number = single(126.87);
binary = float2bin(number);
binary = [repmat('0', 1,31-length(binary)), binary]

%% The bit representation of a floating point number is its own logarithm
% This trick is used to compute the Fast Inverse Square Root (from Quake)
% See https://www.youtube.com/watch?v=p8u_k2LIZyo
mu = 0.0430; % from Quake approximation, but I think this should rather be around 0.0573

log2calculated = log2(number)
log2approx = bin2dec(binary) / 2^23 + mu - 127

exponent_bin = binary(1:8);
mantissa_bin = binary(9:end);
exponent = bin2dec(exponent_bin);
mantissa = bin2dec(mantissa_bin);

manual_float = (1 + mantissa / 2^23) * 2^(exponent - 127)
log2_approx = (mantissa / 2^23) + mu + (exponent - 127) % Quake approximation

%%
x = (0:0.0001:1)';
y = log2(1+x);
plot(x,y, x,x);
axis equal;
xlim([0, 1]);
ylim([0, 1]);

residual = y - x;
mean_residual = mean(residual)
hold on;
plot(x, x+mean_residual);
plot(x, x+0.0430);
hold off;
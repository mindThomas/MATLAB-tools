function bin = float2bin(number)
    num = single(number);
    bin = dec2bin(typecast(num, 'uint32'));
end
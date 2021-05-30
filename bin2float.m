function number = bin2float(bin)    
    num_dec = bin2dec(bin);
    num32 = cast(num_dec, 'uint32');
    number = typecast(num32, 'single');
end
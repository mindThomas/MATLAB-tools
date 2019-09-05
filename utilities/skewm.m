function S = skewm(v)
    if (size(v,1) ~= 3)
        error('skewm: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
    end
    S = zeros(3,3);

    S(1,2) = -v(3,1);
    S(1,3) =  v(2,1);

    S(2,1) =  v(3,1);
    S(2,3) = -v(1,1);

    S(3,1) = -v(2,1);
    S(3,2) =  v(1,1);
end

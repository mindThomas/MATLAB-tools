function v = skewm2vec(S)
    if ( (size(S,1) ~= 3) || (size(S,2) ~= 3) )
        error('skewm2vec: %s', WBM.wbmErrorMsg.WRONG_MAT_DIM);
    end
    v = zeros(3,1);

    v(1,1) = -S(2,3);
    v(2,1) =  S(1,3);
    v(3,1) = -S(1,2);
end

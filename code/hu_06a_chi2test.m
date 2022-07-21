function X2 = hu_06a_chi2test(observed,expected)
    X2 = zeros(width(observed),1);
    
    for i = 1:width(observed)
        for j = 1:height(observed)
            X2(i) = X2(i) + (observed(j,i) - expected(j))^2 / expected(j);
        end
    end
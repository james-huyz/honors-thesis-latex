function odata = hu_06b_gaussamp(orientationdata,i,sigma)
    % Set random seed
    rng('default');
    rng(1);
    
    size = height(orientationdata);
    odata = zeros(size,1);
    x = round(normrnd(i,sigma,round(size*2.25),1));
    x_inrange = x(x >= 0 & x <= 90);
    obls = x_inrange(1:size);
    
    for k = 1:size
        odata(k) = orientationdata(k,obls(k)+1);
    end
end
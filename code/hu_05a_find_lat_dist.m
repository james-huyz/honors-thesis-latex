function lat_dist = hu_05a_find_lat_dist(lats)
    lat_dist = zeros(17,1);
    
    for i = 1:17  % 5*(17-1)+5 = 90°
        freq_at_lat = height(lats(lats > 5*(i-1) & lats <= 5*(i-1) + 5))...
                      / (height(lats)*width(lats));
        lat_dist(i) = freq_at_lat;
    end
end
function mask = myOr(mask1, mask2)
   
    if isempty(mask1)
        if isempty(mask2)
            mask = false;
        else
            mask = mask2;
        end
    elseif isempty(mask2)
        mask = mask1;
    else
        mask = mask1 | mask2;
    end
    
end
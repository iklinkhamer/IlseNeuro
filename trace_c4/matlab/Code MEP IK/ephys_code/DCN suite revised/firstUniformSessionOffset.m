function offset = firstUniformSessionOffset(mouseName)
arguments
    mouseName (1,1) string
end

switch mouseName
    
    case "Dallas"
        offset = 1;  
        
    case "Newark"
        offset = 15;
        
    case "ReserveMouse3"
        offset = 32;
        
    otherwise
        offset = 0;        
        
end
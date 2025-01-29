function type = getMouseType(mname)
    
    mouseList = defaultMice();
    
    mIdx = find(strcmp([mouseList.code], mname), 1); %find(strcmp([mouseList.name], mname), 1);
    
    type = mouseList(mIdx).type;
    
end
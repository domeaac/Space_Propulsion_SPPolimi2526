function BC = initBC(options)
    arguments
        options.T_hotWall = 1800;
        options.T_hotWallFuncHandle = []; 
    end 

    BC.T_hotWall = options.T_hotWall; 
    BC.T_hotWallFuncHandle = options.T_hotWallFuncHandle; 
 


end
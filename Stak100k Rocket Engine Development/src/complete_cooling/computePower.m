function Q = computePower(r_mat,q_vec)

% r_mat 2x(N+1)      first row are the x-coordinates of the nodes of the outer 
%                    boundary, second row are the y-coordinates of the nodes.
% q_vec 1xN          heat flux per unit surface through each element

 N = length(q_vec);

 Q = 0;

 for i = 1 : N 
   
     x1 = r_mat(1, i); x2 = r_mat(1, i+1);
     y1 = r_mat(2, i); y2 = r_mat(2, i+1);
     q_i =q_vec(i);
    
     L_i = sqrt((x2 - x1)^2 + (y2 - y1)^2);
     S_i = pi * (y1 + y2) * L_i; % Superficie tronco di cono
     Q_i =S_i * q_i;

     Q = Q + Q_i;
    
 end

end
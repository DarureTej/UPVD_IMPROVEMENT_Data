  function [distp,distt]=dist_predicition_N(t,N,condition)
% t ->           time instant; 
% N ->           prediciton horizon; 
% condition ->   1 = normalized at operating point; 2 = raw data

    % Loading disturbace data 
    load('dist_solar_lisbon')   

    l    =    length(dist_solar(:,1));
    I    =    dist_solar(t,1);
    Toa  =    dist_solar(t,2);
    Ttin =    dist_solar(t,3);
    
%     condition  = 2;  % [1= with operating point, 2=without operating point]
    switch condition 
        case 1 
            % operating points 
            I_op = 1000;
            Toa_op = 20;
            Ttin_op =10; 

        case 2 
            % without operating points 
            I_op = 1000*0;
            Toa_op = 20*0;
            Ttin_op =10*0; 
    
    end 
    % distt - Disturbance at time t in vector form 
    
    distt = [I-I_op;Toa-Toa_op;Ttin-Ttin_op];

    % distp - Predicted Disturbance for horizon N ...
%                                     at time t in vector form 
    d     =  distt;

    for j=2:N
         d  =  [d  ,[dist_solar((t+j),1)'-I_op;  dist_solar((t+j),2)-Toa_op;   dist_solar((t+j),3)-Ttin_op]];
    end

distp=d;   
  end 


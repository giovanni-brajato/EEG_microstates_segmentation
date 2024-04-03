function [ map ] = topographic_brain_map_64(vector )
%TOPOGRAPHIC_BRAIN_MAP_64 Summary of this function goes here
%   Detailed explanation goes here
map = [0, 0, 0, vector(22), 0, vector(23), 0, vector(24), 0, 0, 0;
                                         0, 0, vector(25:26),0, vector(27), 0, vector(28:29), 0, 0;
                                         0,                       vector(30:38),                          0;
                                         0, vector(39),            vector(1:7),           vector(40), 0 ; 
                                         vector(43), vector(41), vector(8:14), vector(42), vector(44) ;
                                         0, vector(45),          vector(15:21),           vector(46), 0;
                                         0,                         vector(47:55),                         0;
                                         0,0, vector(56:57), 0 , vector(58), 0, vector(59:60), 0, 0;
                                         0, 0,  0,  vector(61), 0, vector(62), 0, vector(63), 0, 0, 0;
                                         0, 0,0 0 0 ,                 vector(64),           0,0,0,0,0;];
%     map = (map + abs(min_potential_value))/(max_potential_value + abs(min_potential_value));                              
%     map = histeq(map);

end


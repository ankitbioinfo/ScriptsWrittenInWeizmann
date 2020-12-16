clear 

sequence={'ZYX', 'ZYZ', 'ZXY', 'ZXZ','YXZ', 'YXY', 'YZX', 'YZY', 'XYZ',...
    'XYX', 'XZY', 'XZX'};


sequence={'XYZ'};

for i=1:length(sequence)
    sequence{i}
dcm = angle2dcm( -pi/3, -pi/4, -pi/6,sequence{i} )
angle=rad2deg(rotm2eul(dcm))
end

R=myroation(-pi/3,-pi/4,-pi/6)

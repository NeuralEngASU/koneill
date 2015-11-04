SimInfo = ParseXML;
MSMSFingerPropAdv(SimInfo,1);

pause(0.1)

prop = [0,1.25,2.5,3.75,5.0];
% prop = [0,NaN,NaN,5,NaN];

MSMSFingerPropAdv(SimInfo,prop);

pause(0.1)

while true
    
    propTmp = prop(end);
    prop(2:end) = prop(1:end-1);
    prop(1) = propTmp;
    
    MSMSFingerPropAdv(SimInfo,prop);
    
    pause(0.1)
    
end
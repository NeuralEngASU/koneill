function [ output_args ] = Demo( m_MSMS )
%DEMO Summary of this function goes here
%   Detailed explanation goes here

m_MSMS = m_MSMS;
end

function tmp = ThumbOut()

for k = 0:0.02:1
   
    tic
    inputStr = ['Finger1=', num2str(k), ';Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;SpherePos1=1.0;SpherePos2=0.0;SpherePos3=0.0;SpherePos4=0.0;SpherePos5=0.0;WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;SphereApp1=0.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

inputStr = 'Finger1=1.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;SpherePos1=0.0;SpherePos2=0.0;SpherePos3=0.0;SpherePos4=0.0;SpherePos5=0.0;WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;';
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);

for k = 1:-0.02:0
   
    tic
    inputStr = ['Finger1=', num2str(k), ';Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;SpherePos1=1.0;SpherePos2=0.0;SpherePos3=0.0;SpherePos4=0.0;SpherePos5=0.0;WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

end

function tmp = IndexOut()

for k = 0:0.02:1
   
    tic
    inputStr = ['Finger1=0.0;Finger2=', num2str(k), ';Finger3=0.0;Finger4=0.0;Finger5=0.0;SpherePos1=0.0;SpherePos2=1.0;SpherePos3=0.0;SpherePos4=0.0;SpherePos5=0.0;WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;SphereApp1=1.0;SphereApp2=0.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

    inputStr = ['Finger1=0.0;Finger2=', num2str(k), ';Finger3=0.0;Finger4=0.0;Finger5=0.0;SpherePos1=0.0;SpherePos2=1.0;SpherePos3=0.0;SpherePos4=0.0;SpherePos5=0.0;WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;SphereApp1=1.0;SphereApp2=0.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);

for k = 1:-0.02:0
   
    tic
    inputStr = ['Finger1=0.0;Finger2=', num2str(k), ';Finger3=0.0;Finger4=0.0;Finger5=0.0;SpherePos1=0.0;SpherePos2=1.0;SpherePos3=0.0;SpherePos4=0.0;SpherePos5=0.0;WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;SphereApp1=1.0;SphereApp2=0.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    

end

function tmp = MiddleOut()

    for k = 0:0.02:1
   
    tic
    inputStr = ['Finger1=0.0;Finger2=', num2str(k), ';Finger3=0.0;Finger4=0.0;Finger5=0.0;SpherePos1=0.0;SpherePos2=1.0;SpherePos3=0.0;SpherePos4=0.0;SpherePos5=0.0;WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;SphereApp1=1.0;SphereApp2=0.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

    inputStr = ['Finger1=0.0;Finger2=', num2str(k), ';Finger3=0.0;Finger4=0.0;Finger5=0.0;SpherePos1=0.0;SpherePos2=1.0;SpherePos3=0.0;SpherePos4=0.0;SpherePos5=0.0;WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;SphereApp1=1.0;SphereApp2=0.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);

for k = 1:-0.02:0
   
    tic
    inputStr = ['Finger1=0.0;Finger2=', num2str(k), ';Finger3=0.0;Finger4=0.0;Finger5=0.0;SpherePos1=0.0;SpherePos2=1.0;SpherePos3=0.0;SpherePos4=0.0;SpherePos5=0.0;WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;SphereApp1=1.0;SphereApp2=0.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

end

end

function tmp = TI()

end

function tmp = TM()

end

function tmp = TIM()

end

% Finger to sphere:      Index, thumb, Index/Middle
% Finger to Thumb Touch: Index/Thumb, Index/Middle
% Have spheres change brightness or opacity when proportions are near.
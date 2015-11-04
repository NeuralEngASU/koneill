function [ ] = Demo( m_MSMS )
%DEMO Summary of this function goes here
%   Detailed explanation goes here

inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=-1.0;SpherePos2=-1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );

pause(3);

IndexOut(m_MSMS);
ThumbOut(m_MSMS);
MiddleOut(m_MSMS);
TIOut(m_MSMS);
TMOut(m_MSMS);
TIMOut(m_MSMS);

m_MSMS.TaskSwitch(1);

TITouch(m_MSMS);
TMTouch(m_MSMS);
TIMTouch(m_MSMS);




end

function [] = ThumbOut(m_MSMS)

pause(0.5);
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=1.0;SpherePos2=-1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=0.5;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);


for k = 0:0.02:1
    
    if k >= 0.9
        brightness = 0;
    else
        brightness = 0.5;
    end
    
    tic
    inputStr = ['Finger1=', num2str(k), ';Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
                'SpherePos1=1.0;SpherePos2=-1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
                'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            	'SphereApp1=', num2str(brightness), ';SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];

    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

pause(0.5);
inputStr = ['Finger1=1.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=0.0;SpherePos2=-1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;'...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=0.5;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);

for k = 1:-0.02:0
    
    if k <= 0.1
        brightness = 0;
    else
        brightness = 0.5;
    end
    
    tic
    inputStr = ['Finger1=', num2str(k), ';Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
                'SpherePos1=0.0;SpherePos2=-1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
                'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            	'SphereApp1=', num2str(brightness), ';SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

pause(0.5)
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=0.0;SpherePos2=-1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );

pause(0.01)
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=-1.0;SpherePos2=-1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );


end

function [] = IndexOut(m_MSMS)

pause(0.5);
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=-1.0;SpherePos2=1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=0.5;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);


for k = 0:0.02:1
    
    if k >= 0.9
        brightness = 0;
    else
        brightness = 0.5;
    end
    
    tic
    inputStr = ['Finger1=0.0;Finger2=', num2str(k), ';Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
                'SpherePos1=-1.0;SpherePos2=1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
                'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
                'SphereApp1=1.0;SphereApp2=', num2str(brightness), ';SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];

    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

pause(0.5);
inputStr = ['Finger1=0.0;Finger2=1.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=-1.0;SpherePos2=0.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=0.5;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);

for k = 1:-0.02:0
    
    if k <= 0.1
        brightness = 0;
    else
        brightness = 0.5;
    end
    
    tic
inputStr = ['Finger1=0.0;Finger2=', num2str(k), ';Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=-1.0;SpherePos2=0.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=', num2str(brightness), ';SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

pause(0.5)
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=-1.0;SpherePos2=0.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );

pause(0.01)
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=-1.0;SpherePos2=-1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );


end

function [] = MiddleOut(m_MSMS)

pause(0.5);
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=-1.0;SpherePos2=-1.0;SpherePos3=1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=0.5;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);


for k = 0:0.02:1
    
    if k >= 0.9
        brightness = 0;
    else
        brightness = 0.5;
    end
    
    tic
    inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=', num2str(k), ';Finger4=0.0;Finger5=0.0;', ...
                'SpherePos1=-1.0;SpherePos2=-1.0;SpherePos3=1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
                'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
                'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=', num2str(brightness), ';SphereApp4=1.0;SphereApp5=1.0;'];

    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

pause(0.5);
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=1.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=-1.0;SpherePos2=-1.0;SpherePos3=0.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1;SphereApp3=0.5;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);

for k = 1:-0.02:0
    
    if k <= 0.1
        brightness = 0;
    else
        brightness = 0.5;
    end
    
    tic
    inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=', num2str(k), ';Finger4=0.0;Finger5=0.0;', ...
                'SpherePos1=-1.0;SpherePos2=-1.0;SpherePos3=0.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
                'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
                'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=', num2str(brightness), ';SphereApp4=1.0;SphereApp5=1.0;'];

    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

pause(0.5)
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=-1.0;SpherePos2=-1.0;SpherePos3=0.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );

pause(0.01)
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=-1.0;SpherePos2=-1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );


end

function [] = TIOut(m_MSMS)

pause(0.5);
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=1.0;SpherePos2=1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=0.5;SphereApp2=0.5;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);


for k = 0:0.02:1
    
    if k >= 0.9
        brightness = 0;
    else
        brightness = 0.5;
    end
    
    tic
    inputStr = ['Finger1=', num2str(k), ';Finger2=', num2str(k), ';Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
                'SpherePos1=1.0;SpherePos2=1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
                'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            	'SphereApp1=', num2str(brightness), ';SphereApp2=', num2str(brightness), ';SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];

    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

pause(0.5);
inputStr = ['Finger1=1.0;Finger2=1.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=0.0;SpherePos2=0.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=0.5;SphereApp2=0.5;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);

for k = 1:-0.02:0
    
    if k <= 0.1
        brightness = 0;
    else
        brightness = 0.5;
    end
    
    tic
    inputStr = ['Finger1=', num2str(k), ';Finger2=', num2str(k), ';Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
                'SpherePos1=0.0;SpherePos2=0.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
                'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            	'SphereApp1=', num2str(brightness), ';SphereApp2=', num2str(brightness), ';SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];

    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

pause(0.5)
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=0.0;SpherePos2=0.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );

pause(0.01)
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=-1.0;SpherePos2=-1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );


end

function [] = TMOut(m_MSMS)

pause(0.5);
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=1.0;SpherePos2=-1.0;SpherePos3=1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=0.5;SphereApp2=1.0;SphereApp3=0.5;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);


for k = 0:0.02:1
    
    if k >= 0.9
        brightness = 0;
    else
        brightness = 0.5;
    end
    
    tic
    inputStr = ['Finger1=', num2str(k), ';Finger2=0.0;Finger3=', num2str(k), ';Finger4=0.0;Finger5=0.0;', ...
                'SpherePos1=1.0;SpherePos2=-1.0;SpherePos3=1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
                'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            	'SphereApp1=', num2str(brightness), ';SphereApp2=1.0;SphereApp3=', num2str(brightness), ';SphereApp4=1.0;SphereApp5=1.0;'];

    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

pause(0.5);
inputStr = ['Finger1=1.0;Finger2=0.0;Finger3=1.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=0.0;SpherePos2=-1.0;SpherePos3=0.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=0.5;SphereApp2=1.0;SphereApp3=0.5;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);

for k = 1:-0.02:0
    
    if k <= 0.1
        brightness = 0;
    else
        brightness = 0.5;
    end
    
    tic
    inputStr = ['Finger1=', num2str(k), ';Finger2=0.0;Finger3=', num2str(k), ';Finger4=0.0;Finger5=0.0;', ...
                'SpherePos1=0.0;SpherePos2=-1.0;SpherePos3=0.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
                'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            	'SphereApp1=', num2str(brightness), ';SphereApp2=1.0;SphereApp3=', num2str(brightness), ';SphereApp4=1.0;SphereApp5=1.0;'];
    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

pause(0.5)
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=0.0;SpherePos2=-1.0;SpherePos3=0.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );

pause(0.01)
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=-1.0;SpherePos2=-1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );


end

function [] = IMOut(m_MSMS)

pause(0.5);
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=-1.0;SpherePos2=1.0;SpherePos3=1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=0.5;SphereApp3=0.5;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);


for k = 0:0.02:1
    
    if k >= 0.9
        brightness = 0;
    else
        brightness = 0.5;
    end
    
    tic
    inputStr = ['Finger1=0.0;Finger2=', num2str(k), ';Finger3=', num2str(k), ';Finger4=0.0;Finger5=0.0;', ...
                'SpherePos1=0.0;SpherePos2=1.0;SpherePos3=1.0;SpherePos4=0.0;SpherePos5=0.0;', ...
                'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
                'SphereApp1=1.0;SphereApp2=', num2str(brightness), ';SphereApp3=', num2str(brightness), ';SphereApp4=1.0;SphereApp5=1.0;'];

    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

pause(0.5);
inputStr = ['Finger1=0.0;Finger2=1.0;Finger3=1.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=-1.0;SpherePos2=0.0;SpherePos3=0.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=0.5;SphereApp3=0.5;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);

for k = 1:-0.02:0
    
    if k <= 0.1
        brightness = 0;
    else
        brightness = 0.5;
    end
    
    tic
    inputStr = ['Finger1=0.0;Finger2=', num2str(k), ';Finger3=', num2str(k), ';Finger4=0.0;Finger5=0.0;', ...
                'SpherePos1=0.0;SpherePos2=0.0;SpherePos3=0.0;SpherePos4=0.0;SpherePos5=0.0;', ...
                'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
                'SphereApp1=1.0;SphereApp2=', num2str(brightness), ';SphereApp3=', num2str(brightness), ';SphereApp4=1.0;SphereApp5=1.0;'];

    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

pause(0.5)
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=-1.0;SpherePos2=-1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );


end

function [] = TIMOut(m_MSMS)

pause(0.5);
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=1.0;SpherePos2=1.0;SpherePos3=1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=0.5;SphereApp2=0.5;SphereApp3=0.5;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);


for k = 0:0.02:1
    
    if k >= 0.9
        brightness = 0;
    else
        brightness = 0.5;
    end
    
    tic
    inputStr = ['Finger1=', num2str(k), ';Finger2=', num2str(k), ';Finger3=', num2str(k), ';Finger4=0.0;Finger5=0.0;', ...
                'SpherePos1=1.0;SpherePos2=1.0;SpherePos3=1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
                'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
                'SphereApp1=', num2str(brightness), ';SphereApp2=', num2str(brightness), ';SphereApp3=', num2str(brightness), ';SphereApp4=1.0;SphereApp5=1.0;'];

    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

pause(0.5);
inputStr = ['Finger1=1.0;Finger2=1.0;Finger3=1.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=0.0;SpherePos2=0.0;SpherePos3=0.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=0.5;SphereApp2=0.5;SphereApp3=0.5;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);

for k = 1:-0.02:0
    
    if k <= 0.1
        brightness = 0;
    else
        brightness = 0.5;
    end
    
    tic
    inputStr = ['Finger1=', num2str(k), ';Finger2=', num2str(k), ';Finger3=', num2str(k), ';Finger4=0.0;Finger5=0.0;', ...
                'SpherePos1=0.0;SpherePos2=0.0;SpherePos3=0.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
                'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
                'SphereApp1=', num2str(brightness), ';SphereApp2=', num2str(brightness), ';SphereApp3=', num2str(brightness), ';SphereApp4=1.0;SphereApp5=1.0;'];

    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

pause(0.5)
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=0.0;SpherePos2=0.0;SpherePos3=0.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );

pause(0.01)
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=-1.0;SpherePos2=-1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );


end

function [] = TITouch(m_MSMS)

pause(0.5);
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=0.6;SpherePos2=1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=0.5;SphereApp2=0.5;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);


for k = 0:0.02:1
    
    if k >= 0.90
        brightness = 0;
    else
        brightness = 0.5;
    end
    
    tic
    inputStr = ['Finger1=', num2str(k*0.6), ';Finger2=', num2str(k), ';Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
                'SpherePos1=0.6;SpherePos2=1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
                'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            	'SphereApp1=', num2str(brightness), ';SphereApp2=', num2str(brightness), ';SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];

    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

pause(0.5);
inputStr = ['Finger1=0.6;Finger2=1.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=0.0;SpherePos2=0.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=0.5;SphereApp2=0.5;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);

for k = 1:-0.02:0
    
    if k <= 0.1
        brightness = 0;
    else
        brightness = 0.5;
    end
    
    tic
    inputStr = ['Finger1=', num2str(k*0.6), ';Finger2=', num2str(k), ';Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
                'SpherePos1=0.0;SpherePos2=0.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
                'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            	'SphereApp1=', num2str(brightness), ';SphereApp2=', num2str(brightness), ';SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];

    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

pause(0.5)
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=0.0;SpherePos2=0.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );

pause(0.01)
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=-1.0;SpherePos2=-1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );


end

function [] = TMTouch(m_MSMS)

pause(0.5);
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=1.0;SpherePos2=-1.0;SpherePos3=1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=0.5;SphereApp2=1.0;SphereApp3=0.5;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);


for k = 0:0.02:1
    
    if k >= 0.90
        brightness = 0;
    else
        brightness = 0.5;
    end
    
    tic
    inputStr = ['Finger1=', num2str(k), ';Finger2=0.0;Finger3=', num2str(k), ';Finger4=0.0;Finger5=0.0;', ...
                'SpherePos1=1.0;SpherePos2=-1.0;SpherePos3=1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
                'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            	'SphereApp1=', num2str(brightness), ';SphereApp2=1.0;SphereApp3=', num2str(brightness), ';SphereApp4=1.0;SphereApp5=1.0;'];

    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

pause(0.5);
inputStr = ['Finger1=1.0;Finger2=0.0;Finger3=1.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=0.0;SpherePos2=-1.0;SpherePos3=0.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=0.5;SphereApp2=1.0;SphereApp3=0.5;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);

for k = 1:-0.02:0
    
    if k <= 0.1
        brightness = 0;
    else
        brightness = 0.5;
    end
    
    tic
    inputStr = ['Finger1=', num2str(k), ';Finger2=0.0;Finger3=', num2str(k), ';Finger4=0.0;Finger5=0.0;', ...
                'SpherePos1=0.0;SpherePos2=-1.0;SpherePos3=0.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
                'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            	'SphereApp1=', num2str(brightness), ';SphereApp2=1.0;SphereApp3=', num2str(brightness), ';SphereApp4=1.0;SphereApp5=1.0;'];
    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

pause(0.5)
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=0.0;SpherePos2=-1.0;SpherePos3=0.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );

pause(0.01)
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=-1.0;SpherePos2=-1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );


end

function [] = TIMTouch(m_MSMS)

pause(0.5);
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=0.75;SpherePos2=1.05;SpherePos3=1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=0.5;SphereApp2=0.5;SphereApp3=0.5;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);


for k = 0:0.02:1
    
    if k >= 0.90
        brightness = 0;
    else
        brightness = 0.5;
    end
    
    tic
    inputStr = ['Finger1=', num2str(k*3/4), ';Finger2=', num2str(k*1.05), ';Finger3=', num2str(k), ';Finger4=0.0;Finger5=0.0;', ...
                'SpherePos1=0.75;SpherePos2=1.05;SpherePos3=1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
                'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
                'SphereApp1=', num2str(brightness), ';SphereApp2=', num2str(brightness), ';SphereApp3=', num2str(brightness), ';SphereApp4=1.0;SphereApp5=1.0;'];

    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

pause(0.5);
inputStr = ['Finger1=0.75;Finger2=1.05;Finger3=1.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=0.0;SpherePos2=0.0;SpherePos3=0.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=0.5;SphereApp2=0.5;SphereApp3=0.5;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );
pause(0.5);

for k = 1:-0.02:0
    
    if k <= 0.1
        brightness = 0;
    else
        brightness = 0.5;
    end
    
    tic
    inputStr = ['Finger1=', num2str(k*3/4), ';Finger2=', num2str(k*1.05), ';Finger3=', num2str(k), ';Finger4=0.0;Finger5=0.0;', ...
                'SpherePos1=0.0;SpherePos2=0.0;SpherePos3=0.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
                'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
                'SphereApp1=', num2str(brightness), ';SphereApp2=', num2str(brightness), ';SphereApp3=', num2str(brightness), ';SphereApp4=1.0;SphereApp5=1.0;'];

    m_MSMS.MoveFingerProp( inputStr );
    timeElap = toc;
    
    if timeElap < 0.01
        pause(0.01 - timeElap)
    end
    
end

pause(0.5)
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=0.0;SpherePos2=0.0;SpherePos3=0.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );

pause(0.01)
inputStr = ['Finger1=0.0;Finger2=0.0;Finger3=0.0;Finger4=0.0;Finger5=0.0;', ...
            'SpherePos1=-1.0;SpherePos2=-1.0;SpherePos3=-1.0;SpherePos4=-1.0;SpherePos5=-1.0;', ...
            'WristRoll=0.5;WristPitch=0.5;WristYaw=0.5;',...
            'SphereApp1=1.0;SphereApp2=1.0;SphereApp3=1.0;SphereApp4=1.0;SphereApp5=1.0;'];
m_MSMS.MoveFingerProp( inputStr );


end






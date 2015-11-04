function [ data ] = GenData(  )

data=ones(250, 500);
for i=1:250
    t = 1:100;
    data(i,50:149) = exp(t/(30+rand(1)*2));

    
    t = 1:300;
    data(i,150:449) = ones(1,300)./exp(t/(90+rand(1)*2))*data(i,149);
    
    
    
end % END FOR
end % END FUNCTION

% EOF
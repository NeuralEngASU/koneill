%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	c2e  (Channel to Electrode)
%		Kevin O'Neill
%		Mario Capecchi Lab
%		20131102
%		v0.1
%		PI: Naveen Nagarajan
%		
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [grid, elec] = c2e( option, arrayType, channels )

%	if nargin < 3
%		fprintf('Not enough input arguments to c2e. Quitting program. DEBUG_nargin\n');
%		return
%	end %END IF
	
%	if option == 1 % (c2e)
%		
%		
%		
%		
%	else if option == 2 % (e2c)
%		
%		
%		
%	end % END IF

	grid = [ 1:10;...
			11:20;...
			21:30;...
			31:40;...
			41:50;...
			51:60;...
			61:70;...
			71:80;...
			81:90;...
			NaN, 91:98, NaN];
			
	elec = [];
	
	
end % END FUNCTION
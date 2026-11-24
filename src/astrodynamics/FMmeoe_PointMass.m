classdef FMmeoe_PointMass < handle % modified Equinoctial orbital elements
	properties
		mu
	end
	methods
		%% FMmeoe_PointMass: function description
		function FMmeoe_PointMass = FMmeoe_PointMass(mu)
			FMmeoe_PointMass.mu = mu;
		end
		
		%% cmpt_KepTerm: function description
		function [dXdt] = cmpt_KepTerm(self, t, x)
            p = x(1);
            f = x(2);
            g = x(3);
            L = x(6); % True longitude
        
            % Auxiliary variables
            q = 1+f*cos(L)+g*sin(L);
        
            dXdt = [	zeros(5,1); ...
						+sqrt(self.mu*p)*((q/p)^2)];
        
		end

		%% cmpt_accel: function description
		function accel = cmpt_accel(self, t, x)
			accel = 0;
		end

		%% cmpt_accel_L: function description
		function accel = cmpt_accel_L(self, t, x)
			accel = self.cmpt_accel(t,x);
		end
		%% cmpt_KepTerm_L: function description
		function KepTerm = cmpt_KepTerm_L(self, t, x)
			KepTerm = self.cmpt_KepTerm(t,x);
		end

	end
end

classdef FMmeoe_OblateEarth_E < handle % modified Equinoctial orbital elements
	properties
        mu
        Re
        J2
	end
	methods
		%% FMmeoe_OblateEarth_E: function description
		function FMmeoe_OblateEarth_E = FMmeoe_OblateEarth_E(mu, Re, J2)
			FMmeoe_OblateEarth_E.mu = mu;
            FMmeoe_OblateEarth_E.Re = Re;
            FMmeoe_OblateEarth_E.J2 = J2;

		end
		
		%% cmpt_KepTerm: function description
		function [dXdt] = cmpt_KepTerm(self, t, x)
            dXdt = zeros(6,1);

		end

		%% cmpt_accel: function description
		function accel = cmpt_accel(self, t, x)
            mEOE_E = x;
            p = mEOE_E(1);
            f = mEOE_E(2);
            g = mEOE_E(3);
            h = mEOE_E(4);
            k = mEOE_E(5);
            L = mEOE_E(6);
            
            q = 1+f*cos(L)+g*sin(L);
            r = p/q;

            b  = 3*self.mu*self.J2*self.Re^2/r^4;
            d = h*sin(L) - k*cos(L);
            q2 = (1 + h^2 + k^2)^2;
            
            aj2r = -1/2*b*( 1 - 12*d^2/q2 );
            aj2t =   -4*b*( d*(h*cos(L) + k*sin(L))/q2 );
            aj2n =   -2*b*( (1 - h^2 - k^2)*d/q2 );

            % Output in Local frame (L)
			accel = [aj2r; aj2t; aj2n];

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

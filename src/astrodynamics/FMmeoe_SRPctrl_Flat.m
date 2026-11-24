classdef FMmeoe_SRPctrl_Flat < handle % % flat-plate model SRP in the Asteroid-centered Hill frame
	properties

		prms
        f_cmptSCattitude
        f_QLS
        checkForEclipse
        Re
        au
    end

	methods
		%% FMmeoe_SRPctrl_Flat: constructor
		function FMmeoe_SRPctrl_Flat = FMmeoe_SRPctrl_Flat(prms, f_cmptSCattitude, f_QLS, checkForEclipse)
            FMmeoe_SRPctrl_Flat.prms 	            = prms;
			FMmeoe_SRPctrl_Flat.f_cmptSCattitude 	= f_cmptSCattitude;
            FMmeoe_SRPctrl_Flat.f_QLS 	            = f_QLS;
			if exist('checkForEclipse', 'var')
				FMmeoe_SRPctrl_Flat.checkForEclipse = checkForEclipse;
                FMmeoe_SRPctrl_Flat.Re = getConstants().Re;
                FMmeoe_SRPctrl_Flat.au = getConstants().planets.e.sma;
			else
				FMmeoe_SRPctrl_Flat.checkForEclipse = false;
                FMmeoe_SRPctrl_Flat.Re = [];
                FMmeoe_SRPctrl_Flat.au = [];
			end

		end
		%% cmpt_KepTerm: function description
		function KepTerm = cmpt_KepTerm(self, t, x)
			KepTerm = [zeros(6,1)];
		end

		%% cmpt_accel: function description
		function [accel, auxout] = cmpt_accel(self, t, x, rkFlag)
            arguments
                self 
                t 
                x 
                rkFlag = true
            end
            % Rotation matrix from Sun (S) to Local frame (L)
            COE = mee2coe(x')';
            QLS = self.f_QLS(t, COE);

            % Check if we're in eclipse condition
            if self.checkForEclipse && self.cmpt_eclipse(x, QLS)
                % In eclipse set everything to zero (ideal coasting):
                % All outputs from attitude function set to zero
                [att, qlawout] = no_qlaw;
                % No solar sail acceleration
                ssacc_L = zeros(3,1);

                % But the mean filter has to be updated as well!
                if self.prms.meanf
                    if rkFlag
                        [~] = rollingMeanRk(t, COE, self.prms.mu);
                    else
                        [~] = rollingMean(t, COE, self.prms.mu);
                    end
                end
                
                % Auxiliary output
                auxout = qlawout;
                auxout.att = att;
                auxout.as_L = ssacc_L;
                auxout.eclipse = true; % Set eclipse variable to true

            else
                % Compute attitude (i.e., cone and clock angles)
			    [att, qlawout] = self.f_cmptSCattitude(t, x, rkFlag);
    
                % Solar-sail acceleration:
                % Acceleration in Sun frame (S)
			    ssacc_S = self.prms.g0 * work_aSRP_hat(att, self.prms.Cprms);
                % Acceleration in Local frame (L)
			    ssacc_L = QLS*ssacc_S;
    
                % Auxiliary output
                auxout = qlawout;
                auxout.att = att;
                auxout.as_L = ssacc_L;
                auxout.eclipse = false; % Set eclipse variable to false

            end

            % output
            accel = ssacc_L;

		end
		%% cmpt_ssacc_L: function description
		function accel = cmpt_accel_L(self, t, x)
			accel = self.cmpt_accel(t,x);
		end
		%% cmpt_KepTerm_L: function description
		function KepTerm = cmpt_KepTerm_L(self, t, x)
			KepTerm = self.cmpt_KepTerm(t,x);
		end

		%% cmpt_eclipse: checks for eclipse using simple cylindrical model
		function eclipse = cmpt_eclipse(self, x, QLS)
            mEOE = x;

            % Sun position unit vector in Local (L) frame
		    rSunHat_L = -QLS(:,1);

            % Eclipse model angles
            r = mee2r(mEOE');
            angEcli = acos(rSunHat_L(1));       % Angle between Earth-sail and Earth-Sun vector
            angEcliSun = acos(self.Re/self.au); % Angle between Earth and Sun
            angEcliSail = acos(self.Re/r);      % Angle between Earth and Earth-Sail vector
   
            % Check angles condition for the cylindrical model
            if (angEcliSun + angEcliSail) < angEcli
			    eclipse = true;
            else
                eclipse = false;
            end

		end



    end

end

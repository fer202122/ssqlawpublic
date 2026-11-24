classdef DynamicsOrbElem < matlab.mixin.Copyable
% classdef DynamicsOrbElem < handle
	properties % Declare class members and methods
		OEtype      % orbital element type
		ForceModels % force models
		bodyModels  % NOT USED
		Nx          % number of elements in state vector
		nFM         % number of force models
		GPE_noGM    % function handle to non-Keplerian dynamics matrix
		GPE         % same as GPE_noGM
		% ImpGPE_noGM
		% ImpGPE
		mu          % gravitational parameter
		% cmpt_dBudx_noGM
		% cmpt_dBudx
		% cmpt_gradGPE_noGM
		% cmpt_HessianGPE_noGM
		% cmpt_rvvec
		% cmpt_drvvec_dxvec
		% cmpt_dxvec_drvvec
		% cmpt_xFromCart
		% cmpt_oe2coe
		% cmpt_coe2oe
		cmpt_dLdt_from_dXdt % sundman transformation factor
		% Rot_I2LVLH
		% Rot_LVLH2I
		eventPrm
		scaleFactor % structure that contains the scaling in longitude 'l' and time 't'
		simOut % simulation output structure
		singularity % NOT USED
		% errorSetting
	end
	methods
		%% DynamicsOrbElem: function description. Constructor
		function DynamicsOrbElem = DynamicsOrbElem(mu, ForceModels, OEtype, scaling, bodyModels, singularity)
			DynamicsOrbElem.mu = mu;
			DynamicsOrbElem.OEtype = OEtype;
			DynamicsOrbElem.nFM = numel(ForceModels);
			if strcmp(ForceModels, 'N/A')
				DynamicsOrbElem.ForceModels = 'N/A';
			else
				DynamicsOrbElem.ForceModels = ForceModels;
			end
			if exist('bodyModels', 'var')
				DynamicsOrbElem.bodyModels = bodyModels;
			else
				DynamicsOrbElem.bodyModels = [];
			end
			if exist('singularity', 'var')
				DynamicsOrbElem.singularity.check = singularity.check;
				if DynamicsOrbElem.singularity.check
					DynamicsOrbElem.singularity.fun = singularity.fun;
				end
			else
				DynamicsOrbElem.singularity.check = 0;
			end


			switch OEtype
                %{
				case 'COE' % classical orbital elements, fast variable: true anomaly.
					DynamicsOrbElem.Nx 					= 6;
					DynamicsOrbElem.GPE_noGM 			= @(GM, x) GaussPlanetaryEqs1(GM, x);
					DynamicsOrbElem.cmpt_dBudx_noGM 	= @(GM, x, u) gradGPE_COE_LVLH(GM, x, u);
					DynamicsOrbElem.cmpt_rvvec 			= @(OE, GM) elem2rv(OE(1,:), OE(2,:), OE(3,:), OE(5,:), OE(4,:), OE(6,:), [], [], [], GM);
					DynamicsOrbElem.cmpt_oe2coe 		= @(t, OE, GM) OE;
					DynamicsOrbElem.Rot_I2LVLH 			= @(OE) Rot_I2LVLH_COE(OE);
					DynamicsOrbElem.Rot_LVLH2I 			= @(OE) Rot_LVLH2I_COE(OE);
					if exist('scaling', 'var') && ~isempty(scaling)
						DynamicsOrbElem.scaleFactor.oe 	= [scaling.l; ones(5,1)];
						DynamicsOrbElem.scaleFactor.l 	= scaling.l;
						DynamicsOrbElem.scaleFactor.t 	= scaling.t;
					else
						DynamicsOrbElem.scaleFactor.oe 	= ones(DynamicsOrbElem.Nx,1);
						DynamicsOrbElem.scaleFactor.l 	= 1;
						DynamicsOrbElem.scaleFactor.t 	= 1;
					end
					DynamicsOrbElem.scaleFactor.L 		= 1;
					DynamicsOrbElem.cmpt_dLdt_from_dXdt = @(X, dXdt) sum(dXdt(4:6));

				case 'mCOE' % classical orbital elements, fast variable: true anomaly. the first variable is p=a(1-e^2) instead of sma
					DynamicsOrbElem.Nx 					= 6;
					DynamicsOrbElem.GPE_noGM 			= @(GM, x) GPE_mCOE_LVLH(GM, x);
					% DynamicsOrbElem.cmpt_dBudx_noGM 	= @(GM, x, u) gradGPE_mCOE_LVLH(GM, x, u);
					DynamicsOrbElem.cmpt_rvvec 			= @(OE, GM) elem2rv(OE(1,:)./(1-OE(2,:).^2), OE(2,:), OE(3,:), OE(5,:), OE(4,:), OE(6,:), [], [], [], GM);
					DynamicsOrbElem.cmpt_oe2coe 		= @(t, OE, GM) [OE(1,:)./(1-OE(2,:).^2); OE(2:6,:)];
					if exist('scaling', 'var') && ~isempty(scaling)
						DynamicsOrbElem.scaleFactor.oe 	= [scaling.l; ones(5,1)];
						DynamicsOrbElem.scaleFactor.l 	= scaling.l;
						DynamicsOrbElem.scaleFactor.t 	= scaling.t;
					else
						DynamicsOrbElem.scaleFactor.oe 	= ones(DynamicsOrbElem.Nx,1);
						DynamicsOrbElem.scaleFactor.l 	= 1;
						DynamicsOrbElem.scaleFactor.t 	= 1;
					end
					DynamicsOrbElem.scaleFactor.L 		= 1;
					DynamicsOrbElem.cmpt_dLdt_from_dXdt = @(X, dXdt) sum(dXdt(4:6));
                %}
				case 'mEOE' % modified equinoctial orbital elements, fast variable: true longitude
					DynamicsOrbElem.Nx 					= 6;
					DynamicsOrbElem.GPE_noGM 			= @(GM, x) GPE_ModEquLVLH(GM, x);
					% DynamicsOrbElem.cmpt_dBudx_noGM 	= @(GM, x, u) cmpt_dBudx_mEOE(GM, x, u);
					% DynamicsOrbElem.cmpt_gradGPE_noGM 	= @(GM, x) gradGPE_ModEquLVLH(GM, x);
					% DynamicsOrbElem.cmpt_HessianGPE_noGM= @(GM, x) second_dGPE_mEOE_LVLH_dMEOE(GM, x);
					% DynamicsOrbElem.cmpt_rvvec 			= @(OE, GM) mEOE2rv(OE(1,:), OE(2,:), OE(3,:), OE(4,:), OE(5,:), OE(6,:), GM);
					% DynamicsOrbElem.cmpt_coe2oe 		= @(COEvec, GM) COE2mEOEvec(COEvec);
					% DynamicsOrbElem.cmpt_xFromCart 		= @(rv, GM) rv2mEOEvec(rv(1:3,:), rv(4:6,:), GM);
					% DynamicsOrbElem.cmpt_drvvec_dxvec 	= @(x, GM) grad_dCart_dmEOE(GM, x);
					% DynamicsOrbElem.cmpt_dxvec_drvvec 	= @(rvvec, GM) grad_dmEOE_dCart(GM, rvvec(1:3,:), rvvec(4:6,:));
					% DynamicsOrbElem.Rot_I2LVLH 			= @(OE) Rot_I2LVLH_MEOE(OE);
					% DynamicsOrbElem.Rot_LVLH2I 			= @(OE) Rot_LVLH2I_MEOE(OE);
                    % Set scaling
					if exist('scaling', 'var') && ~isempty(scaling)
						DynamicsOrbElem.scaleFactor.oe 	= [scaling.l; ones(5,1)];
						DynamicsOrbElem.scaleFactor.l 	= scaling.l;
						DynamicsOrbElem.scaleFactor.t 	= scaling.t;
                        DynamicsOrbElem.scaleFactor.L   = scaling.L;
					else
						DynamicsOrbElem.scaleFactor.oe 	= ones(DynamicsOrbElem.Nx,1);
						DynamicsOrbElem.scaleFactor.l 	= 1;
						DynamicsOrbElem.scaleFactor.t 	= 1;
                        DynamicsOrbElem.scaleFactor.L   = 1;
                    end
                    % Sundman transformation factor
					DynamicsOrbElem.cmpt_dLdt_from_dXdt = @(X, dXdt) dXdt(6);
                %{
				case {'Milan', 'Milan_L'}
					switch OEtype
						case 'Milan' 	% Milankovitch orbital elements, fast variable: mean longitude
						trueOrMeanLon 						= 'meanLon';
						DynamicsOrbElem.ImpGPE_noGM 		= @(GM, x) ImpGPE_Milan(GM, x, trueOrMeanLon);
						DynamicsOrbElem.ImpGPE 				= @(x) DynamicsOrbElem.ImpGPE_noGM(mu, x, trueOrMeanLon);
						case 'Milan_L' 	% Milankovitch orbital elements, fast variable: true longitude
						trueOrMeanLon 						= 'trueLon';
						DynamicsOrbElem.ImpGPE_noGM 		= @(GM, x, dv) ImpGPE_Milan_L(GM, x, dv);
						DynamicsOrbElem.ImpGPE 				= @(x, dv) DynamicsOrbElem.ImpGPE_noGM(GM, x, dv);
						DynamicsOrbElem.scaleFactor.L 		= 1;
						DynamicsOrbElem.cmpt_dLdt_from_dXdt = @(X, dXdt) dXdt(7);
						DynamicsOrbElem.cmpt_drvvec_dxvec 	= @(x, GM) grad_dCart_dMilan_L(GM, x);
					end
					DynamicsOrbElem.Nx 					= 7;
					DynamicsOrbElem.GPE_noGM 			= @(GM, x) GPE_Milan(GM, x, [], [], trueOrMeanLon);
					DynamicsOrbElem.cmpt_dBudx_noGM 	= @(GM, x, u) cmpt_dBudx_milan(GM, x, u, trueOrMeanLon);
					DynamicsOrbElem.cmpt_rvvec 			= @(x, GM) Milan2rv(x(1:3,:), x(4:6,:), x(end,:), GM, trueOrMeanLon);
					DynamicsOrbElem.cmpt_xFromCart 		= @(rv, GM) rv2MilanVec(rv(1:3,:), rv(4:6,:), GM, trueOrMeanLon);

					if exist('scaling', 'var') && ~isempty(scaling)
						DynamicsOrbElem.scaleFactor.oe 	= [scaling.l^2/scaling.t * ones(3,1); ones(4,1)];
						DynamicsOrbElem.scaleFactor.l 	= scaling.l;
						DynamicsOrbElem.scaleFactor.t 	= scaling.t;
					else
						DynamicsOrbElem.scaleFactor.oe 	= ones(DynamicsOrbElem.Nx,1);
						DynamicsOrbElem.scaleFactor.l 	= 1;
						DynamicsOrbElem.scaleFactor.t 	= 1;
                    end
                %}
            end
            % General scaling
			% DynamicsOrbElem.scaleFactor.v 				= DynamicsOrbElem.scaleFactor.l / DynamicsOrbElem.scaleFactor.t;
			% DynamicsOrbElem.scaleFactor.a 				= DynamicsOrbElem.scaleFactor.v / DynamicsOrbElem.scaleFactor.t;
			% DynamicsOrbElem.scaleFactor.GM 				= DynamicsOrbElem.scaleFactor.l^3 / DynamicsOrbElem.scaleFactor.t^2;
			DynamicsOrbElem.scaleFactor.dxdt 			= DynamicsOrbElem.scaleFactor.oe / DynamicsOrbElem.scaleFactor.t;
			DynamicsOrbElem.scaleFactor.dxdL 			= DynamicsOrbElem.scaleFactor.oe / DynamicsOrbElem.scaleFactor.L;
			% DynamicsOrbElem.scaleFactor.dxdt_dx 		= DynamicsOrbElem.scaleFactor.dxdt ./ (DynamicsOrbElem.scaleFactor.oe');
			% DynamicsOrbElem.scaleFactor.A_2nd 			= repmat(DynamicsOrbElem.scaleFactor.dxdt_dx, 1, 1, DynamicsOrbElem.Nx) ./ reshape(DynamicsOrbElem.scaleFactor.oe, 1, 1, DynamicsOrbElem.Nx);
			% DynamicsOrbElem.scaleFactor.XwSTM 			= [DynamicsOrbElem.scaleFactor.oe; vec(DynamicsOrbElem.scaleFactor.oe./(DynamicsOrbElem.scaleFactor.oe'))];

			DynamicsOrbElem.GPE 						= @(x) DynamicsOrbElem.GPE_noGM(mu, x); % Computes the non-Keplerian dyamics matrix
			% DynamicsOrbElem.cmpt_dBudx 					= @(x, u) DynamicsOrbElem.cmpt_dBudx_noGM(mu, x, u);
        end
        %{
		%% resetParam: function description
		function resetParam(self, prmname, prmval)
			if isfield(self, prmname)
				self.prmname = prmval;
				self.GPE = @(x) self.GPE_noGM(self.mu, x);
				self.ImpGPE = @(x) self.ImpGPE_noGM(self.mu, x);
				self.cmpt_dBudx = @(x, u) self.cmpt_dBudx_noGM(self.mu, x, u);
			end
			for ii = 1:self.nFM
				if isprop(self.ForceModels{ii}, prmname)
					self.ForceModels{ii}.(prmname) = prmval;
				end
			end
		end
		%% setErrorSettings: function description
		function setErrorSettings(self, prmname, prmval)
			for ii=1:numel(prmname)
				self.errorSetting.(prmname{ii}) = prmval{ii};
			end
		end
        %}
		%% ODE: function description
		function [total_dXdt] = ODE(self, t, X)
			total_accel 	= zeros(3,1);
			total_KepTerm 	= zeros(self.Nx,1);

			t = t * self.scaleFactor.t;
			X = X.* self.scaleFactor.oe;
			% X(end) = wrapTo2Pi(X(end));
			for ii = 1:self.nFM
				total_KepTerm = total_KepTerm + self.ForceModels{ii}.cmpt_KepTerm(t, X);
				total_accel = total_accel + self.ForceModels{ii}.cmpt_accel(t, X);
            end
			total_dXdt = (total_KepTerm + self.GPE(X) * total_accel) ./ self.scaleFactor.dxdt;

		end
        %{
		%% ODE_externalAccel: function description
		function [total_dXdt] = ODE_externalAccel(self, t, X, accel_nd)
			total_accel 	= accel_nd * (self.scaleFactor.v / self.scaleFactor.t);
			total_KepTerm 	= zeros(self.Nx,1);

			t = t * self.scaleFactor.t;
			X = X.* self.scaleFactor.oe;
			for ii = 1:self.nFM
				total_KepTerm = total_KepTerm + self.ForceModels{ii}.cmpt_KepTerm(t, X);
				total_accel = total_accel + self.ForceModels{ii}.cmpt_accel(t, X);
			end
			total_dXdt = (total_KepTerm + self.GPE(X) * total_accel) ./ self.scaleFactor.dxdt;
		end
        %}
		%% ODE_TL: function description
		function [total_dXdL] = ODE_TL(self, L, X)
			% ODE in true longitude
			total_accel 	= zeros(3,1);
			total_KepTerm 	= zeros(self.Nx,1);

			L = L * self.scaleFactor.L;
			t = X(end) * self.scaleFactor.t;
			OE = X(1:self.Nx).* self.scaleFactor.oe;
            OE(end) = L; % Store independent variable in OE so that it reaches the acceleration functions
			for ii = 1:self.nFM
				total_KepTerm = total_KepTerm + self.ForceModels{ii}.cmpt_KepTerm_L(t, OE);
				total_accel = total_accel + self.ForceModels{ii}.cmpt_accel_L(t, OE);
            end
			total_dOEdt = total_KepTerm + self.GPE(OE) * total_accel;
			dtdL 		= 1 / self.cmpt_dLdt_from_dXdt(OE, total_dOEdt); % Sundman transformation factor
            % Assemble total derivative and scale it
			total_dXdL  = [total_dOEdt * dtdL ./ self.scaleFactor.dxdL; dtdL / (self.scaleFactor.t / self.scaleFactor.L)];

		end
        %{
		%% ODEwithSTM: function description
		function [total_dXdt_wSTM] = ODEwithSTM(self, t, X)
			state = X(1:self.Nx);
			STMvec = X((self.Nx+1):(self.Nx+self.Nx^2));

			t = t * self.scaleFactor.t;
			state = state.* self.scaleFactor.oe;

			total_accel 	= zeros(3,1);
			total_KepTerm 	= zeros(self.Nx,1);
			total_dkdx 		= zeros(self.Nx,self.Nx);
			total_dudx 		= zeros(3,self.Nx);
			for ii = 1:self.nFM
				total_KepTerm 	= total_KepTerm + self.ForceModels{ii}.cmpt_KepTerm(t, state);
				total_accel 	= total_accel + self.ForceModels{ii}.cmpt_accel(t, state);
				total_dkdx 		= total_dkdx + self.ForceModels{ii}.cmpt_dkdx(t,state);
				total_dudx 		= total_dudx + self.ForceModels{ii}.cmpt_dudx(t,state);
			end
			Bmat 			= self.GPE(state);
			total_dXdt 		= (total_KepTerm + Bmat * total_accel) ./ self.scaleFactor.dxdt;
			if sum(abs(total_accel))==0
				total_dXdt_dX 	= (total_dkdx + Bmat * total_dudx) ./ self.scaleFactor.dxdt_dx;
			else
				total_dXdt_dX 	= (total_dkdx + self.cmpt_dBudx(state, total_accel) + Bmat * total_dudx) ./ self.scaleFactor.dxdt_dx;
			end

			dSTMdt = total_dXdt_dX * reshape(STMvec, self.Nx, self.Nx);
			total_dXdt_wSTM = [total_dXdt; dSTMdt(:)];
		end
		
		%% ODEwithSTT2nd: function description
		function [total_dXdt_wSTT] = ODEwithSTT2nd(self, t, X)
			state = X(1:self.Nx);

			% reshape STM
			phi_first 	= reshape(X((self.Nx+1):(self.Nx+self.Nx^2)), self.Nx, self.Nx);

			% reshape second order STT
			second_start = self.Nx+self.Nx^2;
			phi_second 	= reshape(X((second_start+1):(second_start+self.Nx^3)), self.Nx, self.Nx, self.Nx);

			% compute necessary terms
			t = t * self.scaleFactor.t;
			state = state.* self.scaleFactor.oe;

			total_accel 		= zeros(3,1);
			total_KepTerm 		= zeros(self.Nx,1);
			total_dKep_1st 		= zeros(self.Nx,self.Nx);
			total_dKep_2nd 		= zeros(self.Nx,self.Nx,self.Nx);
			total_accel_A_1st 	= zeros(3,self.Nx);
			total_accel_A_2nd 	= zeros(3,self.Nx,self.Nx);
			for ii = 1:self.nFM
				% Keplerian terms
				total_KepTerm 		= total_KepTerm + self.ForceModels{ii}.cmpt_KepTerm(t, state);
				total_dKep_1st 		= total_dKep_1st + self.ForceModels{ii}.cmpt_dkdx(t,state);
				total_dKep_2nd 		= total_dKep_2nd + self.ForceModels{ii}.cmpt_dKep_2nd(t,state);

				% acceleration terms
				total_accel 		= total_accel + self.ForceModels{ii}.cmpt_accel(t, state);
				total_accel_A_1st 	= total_accel_A_1st + self.ForceModels{ii}.cmpt_dudx(t,state);
				total_accel_A_2nd 	= total_accel_A_2nd + self.ForceModels{ii}.cmpt_A_2nd(t,state);
			end

			% orbital element variation
			Bmat 			= self.GPE(state);
			total_dXdt 		= (total_KepTerm + Bmat * total_accel) ./ self.scaleFactor.dxdt;


			% first order
			if sum(abs(total_accel))==0
				% if no perturbation, use just partials of Keplerian term
				total_A_first 	=  total_dKep_1st ./ self.scaleFactor.dxdt_dx;
			else
				total_A_first 	= (total_dKep_1st + self.cmpt_dBudx(state, total_accel) + Bmat * total_accel_A_1st) ./ self.scaleFactor.dxdt_dx;
			end
			phiDot_first = total_A_first * phi_first;


			% second order
			if sum(abs(total_accel))==0
				total_A_second = total_dKep_2nd ./ self.scaleFactor.A_2nd;
			else
				t11 = zeros(self.Nx,self.Nx,self.Nx);
				t13 = zeros(self.Nx,self.Nx,self.Nx);
				t33 = zeros(self.Nx,self.Nx,self.Nx);
				% all acceleration terms are expressed in the RTN frame
				grad_GPE = self.cmpt_gradGPE_noGM(self.mu, state);
				for ii = 1:self.Nx
					for g1 = 1:3
						for aa = 1:self.Nx          
							for bb = 1:self.Nx          
								t11(ii,aa,bb) = Bmat(ii,g1) * total_accel_A_2nd(g1,aa,bb);
								t13(ii,aa,bb) = grad_GPE(ii,g1,aa) * total_accel_A_1st(g1,bb);
								t33(ii,aa,bb) = grad_GPE_2nd(ii,g1,aa,bb) * total_accel(g1);
							end
						end
					end
				end
				t31 = t13;
				total_A_second 	= (t11 + t13 + t31 + t33)./ self.scaleFactor.A_2nd;
				% total_A_second 	= t11 + t12 + t13 + t21 + t22 + t23 + t31 + t32 + t33;
			end

			phiDot_second = zeros(self.Nx,self.Nx,self.Nx);
			for ii = 1:self.Nx
				for jj = 1:self.Nx
					for kk = 1:self.Nx          
						for ll = 1:self.Nx
							phiDot_second(ii,jj,kk) = phiDot_second(ii,jj,kk) + total_A_first(ii,ll)*phi_second(ll,jj,kk);
							for mm = 1:self.Nx
								phiDot_second(ii,jj,kk) = phiDot_second(ii,jj,kk) + total_A_second(ii,ll,mm)*phi_first(ll,jj)*phi_first(mm,kk);
							end
						end
					end
				end
			end

			% sum up
			total_dXdt_wSTT = [total_dXdt; phiDot_first(:); phiDot_second(:)];
		end
		

		%% ODE_ImpGPE: function description
		function Delta_X = ODE_ImpGPE(self, X, DV, N_order)
			if N_order == 2
				[B1mat, B2mat] = self.ImpGPE(X);
				quadTerm = [0; 0; 0; ...
							DV'*B2mat(:,:,1)*DV; ...
							DV'*B2mat(:,:,2)*DV; ...
							DV'*B2mat(:,:,3)*DV; ...
							0	];
			elseif N_order == 1
				B1mat = self.ImpGPE(X);
				quadTerm = 0;
			end
			linTerm = B1mat * DV;
			Delta_X = linTerm + quadTerm;
		end

		%% cmptElemAfterStatDV_LinFB: function description
		function [OE_p, dxp_dxm, dxp_dparam] = cmptElemAfterStatDV_LinFB(self, OEm, dv_det, Kgain, OE_mean, dimOrNondim, wrapFastVar)

			if isempty(Kgain)
				% compute TCM -- deterministic
				dv_TCM = dv_det;
			else
				% compute TCM -- linear FB
				dv_TCM = dv_det + Kgain * (OEm - OE_mean);
			end
		
			if nargout==1
				% compute element mapping
				[OE_p] = self.cmptElemAfterDV(OEm, dv_TCM, dimOrNondim, wrapFastVar);

			elseif nargout>1
				% analytical gradients
				[OE_p, dxp_dxm1, dxp_ddv] = self.cmptElemAfterDV(OEm, dv_TCM, dimOrNondim, wrapFastVar);

				if isempty(Kgain)
					% for deterministic TCM
					ddv_ddvdet = eye(3);
					ddv_dprm = ddv_ddvdet;

				else
					% for statistical TCM
					ddv_dxm = Kgain;
					ddv_ddvdet = eye(3);
					ddv_dK = kron((OEm - OE_mean)', eye(3));
					ddv_dxnom = - Kgain;
					ddv_dprm = [ddv_ddvdet, ddv_dK, ddv_dxnom];
				end
			
				% summarize
				dxp_dxm = dxp_dxm1 + dxp_ddv * ddv_dxm;
				dxp_dparam = dxp_ddv * ddv_dprm;
			end
		end

		%% cmptElemAfterDV: function description
		function [OE_p, grad_OE, grad_dv] = cmptElemAfterDV(self, OE_m, dv, dimOrNondim, wrapFastVar)
			if strcmp(self.OEtype, 'Milan_L')
				if dimOrNondim==1 % dimensional
					[D_OE, dgdx, dgdv] = ImpGPE_Milan_L(self.mu, OE_m, dv);
				elseif dimOrNondim==2
					GM_nd 			= self.mu / self.scaleFactor.GM;
					[D_OE, dgdx, dgdv] = ImpGPE_Milan_L(GM_nd, OE_m, dv);
				end
				D_OE(end) = wrapFastVariable(D_OE(end));
				OE_p 	= OE_m + D_OE;
				grad_OE = eye(self.Nx) + dgdx;
				grad_dv = dgdv;
			else
				if dimOrNondim==1 % dimensional
					[rvec, vvec_m] 	= self.cmpt_rvvec(OE_m, self.mu);
					rvvec_p 		= [rvec; vvec_m + dv];
					OE_p 			= self.cmpt_xFromCart(rvvec_p, self.mu);

				elseif dimOrNondim==2
					GM_nd 			= self.mu / self.scaleFactor.GM;
					[rvec, vvec_m] 	= self.cmpt_rvvec(OE_m, GM_nd);
					rvvec_p 		= [rvec; vvec_m + dv];
					OE_p 			= self.cmpt_xFromCart(rvvec_p, GM_nd);
				end
				D_OE = cmpt_OEdifference(OE_p, OE_m);
				OE_p = OE_m + D_OE;

				% compute analytical gradient
				if nargout>1
					if dimOrNondim==1 % dimensional
						drvm_dOEm 	= self.cmpt_drvvec_dxvec(OE_m, self.mu);
						drvp_drvm 	= eye(6);
						dOEp_drvp 	= self.cmpt_dxvec_drvvec(rvvec_p, self.mu);
						drvp_ddv 	= [zeros(3,3); eye(3)];

					elseif dimOrNondim==2
						drvm_dOEm 	= self.cmpt_drvvec_dxvec(OE_m, GM_nd);
						drvp_drvm 	= eye(6);
						dOEp_drvp 	= self.cmpt_dxvec_drvvec(rvvec_p, GM_nd);
						drvp_ddv 	= [zeros(3,3); eye(3)];

					end
					grad_OE = dOEp_drvp * drvp_drvm * drvm_dOEm;
					grad_dv = dOEp_drvp * drvp_ddv;
				end
			end
			if wrapFastVar
				OE_p(end) 		= wrapFastVariable(OE_p(end));
			end
		end

		%% BplaneMapping: function description
		function [output, grad] = BplaneMapping(self, t, OE, bodyID, outType, dimOrNondim)
			BplaneAnalyticGradient = 1;
			if dimOrNondim==1
				sf_Bprm = ones(4,1);
			elseif dimOrNondim==2
				sf_Bprm = [	self.scaleFactor.l; ...
							self.scaleFactor.l; ...
							self.scaleFactor.t; ...
							self.scaleFactor.l ...
							];
			end

			if ~exist('outType', 'var')
				outType = 1;
			end
			if dimOrNondim==1 % dimensional
				t_dim 			= t;
				[rvec, vvec] 	= self.cmpt_rvvec(OE, self.mu);
			elseif dimOrNondim==2
				t_dim 			= t * self.scaleFactor.t;
				GM_nd 			= self.mu / self.scaleFactor.GM;
				[rvec, vvec] 	= self.cmpt_rvvec(OE, GM_nd);
				rvec 			= rvec * self.scaleFactor.l;
				vvec 			= vvec * self.scaleFactor.v;
			end
			GM_planet = self.bodyModels.Sbodies{bodyID}.mu;
			rv_planet = self.bodyModels.get_Xpj(t_dim, bodyID);
			
			r_rel = rvec - rv_planet(1:3);
			v_rel = vvec - rv_planet(4:6);
			
			if BplaneAnalyticGradient
				[BdotT, BdotR, LTOF, RadPeri, Psi, dBprm_drvvec] = cmptBplaneParams(GM_planet, r_rel, v_rel, 3);
			else
				[BdotT, BdotR, LTOF, RadPeri, Psi] = cmptBplaneParams(GM_planet, r_rel, v_rel, 3);
			end
			BdotT 	= BdotT 	/ sf_Bprm(1);
			BdotR 	= BdotR 	/ sf_Bprm(2);
			LTOF 	= LTOF 		/ sf_Bprm(3);
			RadPeri = RadPeri 	/ sf_Bprm(4);
			
			if outType==1
				% output only the B-plane location
				output = [BdotT; BdotR];
			elseif outType==2
				% output the B-plane params: B-plane location + LTOF
				output = [BdotT; BdotR; LTOF];
			elseif outType==0
				output = [BdotT; BdotR; RadPeri; Psi];
			end

			if nargout>1 % analytical gradient
				if dimOrNondim==1 % dimensional
					drvvec_dOE = self.cmpt_drvvec_dxvec(OE, self.mu);
				elseif dimOrNondim==2
					drvvec_dOE = self.cmpt_drvvec_dxvec(OE.*self.scaleFactor.oe, self.mu);
				end

				% numerical gradient: dimensional
				if ~BplaneAnalyticGradient
					dBprm_drvvec = zeros(4,6);
					pertVec = 1E-8 * norm(r_rel,2) * eye(3);
					for ii=1:3
						[BdotT_p, BdotR_p, LTOF_p, RadPeri_p, Psi_p] = cmptBplaneParams(GM_planet, r_rel + pertVec(:,ii), v_rel, 1);
						dBprm = [	BdotT_p 	- BdotT 	* sf_Bprm(1); ...
									BdotR_p 	- BdotR 	* sf_Bprm(2); ...
									LTOF_p 		- LTOF 		* sf_Bprm(3); ...
									RadPeri_p 	- RadPeri 	* sf_Bprm(4)];
						dBprm_drvvec(:,ii) = dBprm/pertVec(ii,ii);
					end
					pertVec = 1E-8 * norm(v_rel,2) * eye(3);
					for ii=4:6
						[BdotT_p, BdotR_p, LTOF_p, RadPeri_p, Psi_p] = cmptBplaneParams(GM_planet, r_rel, v_rel + pertVec(:,ii-3), 1);
						dBprm = [	BdotT_p 	- BdotT 	* sf_Bprm(1); ...
									BdotR_p 	- BdotR 	* sf_Bprm(2); ...
									LTOF_p 		- LTOF 		* sf_Bprm(3); ...
									RadPeri_p 	- RadPeri 	* sf_Bprm(4)];
						dBprm_drvvec(:,ii) = dBprm/pertVec(ii-3,ii-3);
					end
				end

				% compute partial of the planet state w.r.t. time
				dt = 1E-8 * 86400;
				rvp_ = self.bodyModels.get_Xpj(t_dim + dt, bodyID);
				drvp_dt = (rvp_ - rv_planet) / dt;

				% summary: dimensional
				dBT_dOE 	= dBprm_drvvec(1,:) * drvvec_dOE;
				dBR_dOE 	= dBprm_drvvec(2,:) * drvvec_dOE;
				dLTOF_dOE 	= dBprm_drvvec(3,:) * drvvec_dOE;
				dRp_dOE 	= dBprm_drvvec(4,:) * drvvec_dOE;

				dBT_dt 		= -dBprm_drvvec(1,:) * drvp_dt;
				dBR_dt 		= -dBprm_drvvec(2,:) * drvp_dt;
				dLTOF_dt 	= -dBprm_drvvec(3,:) * drvp_dt;
				dRp_dt 		= -dBprm_drvvec(4,:) * drvp_dt;

				if dimOrNondim==1 % dimensional
					grad = struct(	'dBT_dOE', dBT_dOE, ...
									'dBR_dOE', dBR_dOE, ...
									'dLTOF_dOE', dLTOF_dOE, ...
									'dRp_dOE', dRp_dOE, ...
									'dBT_dt', dBT_dt, ...
									'dBR_dt', dBR_dt, ...
									'dLTOF_dt', dLTOF_dt, ...
									'dRp_dt', dRp_dt);
				elseif dimOrNondim==2
					grad = struct(	'dBT_dOE', dBT_dOE ./ (self.scaleFactor.l./(self.scaleFactor.oe')), ...
									'dBR_dOE', dBR_dOE ./ (self.scaleFactor.l./(self.scaleFactor.oe')), ...
									'dLTOF_dOE', dLTOF_dOE ./ (self.scaleFactor.t./(self.scaleFactor.oe')), ...
									'dRp_dOE', dRp_dOE ./ (self.scaleFactor.l./(self.scaleFactor.oe')), ...
									'dBT_dt', dBT_dt ./ (self.scaleFactor.l/(self.scaleFactor.t)), ...
									'dBR_dt', dBR_dt ./ (self.scaleFactor.l/(self.scaleFactor.t)), ...
									'dLTOF_dt', dLTOF_dt ./ (self.scaleFactor.t/(self.scaleFactor.t)), ...
									'dRp_dt', dRp_dt ./ (self.scaleFactor.l/(self.scaleFactor.t)));
				end
			end

		end		

		%% cmptTCM_BplaneTarget: function description
		function [dvTCM, LTOF] = cmptTCM_BplaneTarget(self, t, OE_m, bodyID, B_target, dv_det, alg, dimOrNondim)
			% alg = 1; % target B-plane location + LTOF,	finite differentiation
			% alg = 2; % target B-plane location, 			finite differentiation
	
			sf_v 		= self.scaleFactor.l / self.scaleFactor.t;
			sf_Bvec 	= [self.scaleFactor.l; self.scaleFactor.l; self.scaleFactor.t];
	
			if dimOrNondim==1
				t_dim 			= t;
				[rvec, vvec] 	= self.cmpt_rvvec(OE_m, self.mu);

			elseif dimOrNondim==2
				t_dim 			= t * self.scaleFactor.t;
				GM_nd 			= self.mu / self.scaleFactor.GM;
				[rvec, vvec] 	= self.cmpt_rvvec(OE_m, GM_nd);
				rvec 			= rvec * self.scaleFactor.l;
				vvec 			= vvec * self.scaleFactor.v;

			end
			GM_planet = self.bodyModels.Sbodies{bodyID}.mu;
			rv_planet = self.bodyModels.get_Xpj(t_dim, bodyID);

			r_ref = rvec - rv_planet(1:3);
			v_ref = vvec - rv_planet(4:6) + dv_det;
			
			% compute B-plane param w/out TCM
			[B_T, B_R, LTOF] = cmptBplaneParams(GM_planet, r_ref, v_ref, 1);
			if alg==1
				B_ref = [B_T; B_R; LTOF];
				B_target = B_target(1:3);
			elseif alg==2  % 
				B_ref = [B_T; B_R];
				B_target = B_target(1:2);
				sf_Bvec(end) = [];
			end
			
			% differential corrector
			eps_v = 1E-08 * sf_v;
			dif = 1; tol = 1e-6; ii = 1; iMax = 100;
			dV = zeros(3,iMax);
			
			isterminate = 0;
			while ~isterminate
				[dB_dV] 	= cmptBplanePartial(GM_planet, r_ref, v_ref, B_ref, eps_v, alg);
				[dB_dV_] 	= cmptBplanePartial(GM_planet, r_ref, v_ref, B_ref, eps_v*10, alg);
			
				if norm(dB_dV - dB_dV_)/norm(dB_dV) > 0.01
					error('V_perturbation is not appropriate')
					% eps_v = eps_v/10;
				else
					ii = ii + 1;
			
					B_err = B_target - B_ref;
					delta_V = dB_dV\B_err;
					% delta_V = inv(dB_dV)*B_err;
					if alg==1
						% do nothing
					elseif alg==2
						delta_V = [delta_V; 0];
					end
			
					v_ref = v_ref + delta_V;
			
					[B_T, B_R, LTOF] = cmptBplaneParams(GM_planet, r_ref, v_ref, 1);
					if alg==1
						B_ref = [B_T; B_R; LTOF];
					elseif alg ==2
						B_ref = [B_T; B_R];
					end
				
					dV(:,ii) = delta_V;
					dif = norm(B_err./sf_Bvec,2);
				end
				if dif < tol || iMax < ii
					isterminate = 1;
				end
			end
			delta_V = sum(dV,2);
			
			dvTCM = dv_det + delta_V;

			if dimOrNondim==2
				dvTCM = dvTCM / self.scaleFactor.v;
				LTOF = LTOF / self.scaleFactor.t;
			end
		end
		
        %}
		%% PropagateSC: function description
		function [t_his, X, sol] = PropagateSC(self, T_span, X0, options, odeSolver)
			if ~exist('odeSolver')
				odeSolver = @(ode, tspan, x0, opts) ode45(ode, tspan, x0, opts);
			end
			self.scaleFactor.iv = self.scaleFactor.t;
			self.scaleFactor.y = self.scaleFactor.oe;

			ode = @(t,X) self.ODE(t, X);
			sol = odeSolver(ode, T_span / self.scaleFactor.iv, X0 ./ self.scaleFactor.y, options);
			self.simOut.sol = sol;
			t_his = sol.x;
			X = self.getTrajectory(sol, t_his);
		end
        %{
		%% PropagateSCwithSTM: function description
		function [t_his, XwSTM, sol] = PropagateSCwithSTM(self, T_span, X0, options)
			if ~exist('odeSolver')
				odeSolver = @(ode, tspan, x0, opts) ode45(ode, tspan, x0, opts);
			end
			self.scaleFactor.iv = self.scaleFactor.t;
			self.scaleFactor.y = self.scaleFactor.XwSTM;

			ode = @(t,X) self.ODEwithSTM(t, X);
			sol = odeSolver(ode, T_span / self.scaleFactor.iv, X0 ./ self.scaleFactor.y, options);
			self.simOut.sol = sol;
			t_his = sol.x;
			XwSTM = self.getTrajectoryWSTM(sol, t_his);
		end
        %}
		%% PropagateSC_TL: function description
		function [t_his, oe_his, sol] = PropagateSC_TL(self, Lspan, OE0, options, odeSolver)
			% propagate in true longitude
			if ~exist('odeSolver')
				odeSolver = @(ode, Lspan, x0, opts) ode45(ode, Lspan, x0, opts);
            end
            % Scale factor independent variable
			self.scaleFactor.iv = self.scaleFactor.L;
            % Scale factor state
			self.scaleFactor.y = [self.scaleFactor.oe; self.scaleFactor.t];

			ode = @(t,X) self.ODE_TL(t, X);

			X0 = [OE0; 0]; % orbital elements + time
			sol = odeSolver(ode, Lspan / self.scaleFactor.iv, X0 ./ self.scaleFactor.y, options);
			self.simOut.sol = sol;
			sol.independentVar = 'anomaly';
			TL_his = sol.x;
			[oe_his, t_his] = self.getTrajectory(sol, TL_his);
		end

		%% getTrajectory: function description
		function [oe_his, t_his] = getTrajectory(self, sol, x_his, showMsg)
			if ~exist('showMsg', 'var')
				showMsg = 0;
			end
			Nt = length(x_his);
			msg = ['The sol data format is not compatible with "deval" function. Output with a simple interpolation'];
			if ~isfield(sol, 'independentVar')
				sol.independentVar = 'time';
			end
			switch sol.independentVar
				case 'time'
				if isfield(sol, 'ie') && ~isempty(sol.ie)
					t_his = linspace(x_his(1), sol.xe*0.999, Nt) * self.scaleFactor.t;
				else
					t_his = linspace(x_his(1), x_his(end)*0.999, Nt) * self.scaleFactor.t;
				end
				try
					oe_his = deval(sol, t_his / self.scaleFactor.iv) .* repmat(self.scaleFactor.y, 1, Nt);
				catch
					if showMsg
						disp(msg)
					end
					t_his_ode = sol.x * self.scaleFactor.t;
					oe_his_ode = sol.y .* repmat(self.scaleFactor.y, 1, length(t_his_ode));
					t_his = linspace(t_his_ode(1), t_his_ode(end), Nt);
					oe_his = InterpArrays(t_his_ode', oe_his_ode', t_his, 'interp1')';
				end

				case 'anomaly'
				if sol.ie==1
					L_his = linspace(x_his(1), sol.xe*0.999, Nt) * self.scaleFactor.iv;
				else
					L_his = linspace(x_his(1), x_his(end)*0.999, Nt) * self.scaleFactor.iv;
				end
				try
					y_his = deval(sol, L_his / self.scaleFactor.iv) .* repmat(self.scaleFactor.y, 1, Nt);
					t_his = y_his(end,:);
					oe_his = y_his(1:self.Nx,:);
				catch
					if showMsg
						disp(msg)
					end
					L_his_sol = sol.x * self.scaleFactor.iv;
					y_his_sol = sol.y .* repmat(self.scaleFactor.y, 1, length(L_his_sol));

					isIncreasing = all(diff(sol.x)>0);
					if ~isIncreasing
						% reverse the order for interp1
						L_his_sol = flip(L_his_sol);
						y_his_sol = flip(y_his_sol,2);
					end
					y_his = InterpArrays(L_his_sol', y_his_sol', L_his, 'interp1')';
					t_his = y_his(end,:);

					if ~isIncreasing
						% reverse again the order to make it in the original direction
						t_his = flip(t_his);
						y_his = flip(y_his,2);
					end

					% idxnan = find(isnan(t_his));
					% t_his(idxnan) = [];
					% y_his(:,idxnan) = [];
					oe_his = y_his(1:self.Nx,:);
				end

			end
		end		
        %{
		%% getTrajectoryWSTM: function description
		function [XwSTM_his, t_his] = getTrajectoryWSTM(self, sol, t_his)
			if sol.ie==1
				t_his = linspace(t_his(1), sol.xe*0.999, length(t_his)) * self.scaleFactor.iv;
			end
			XwSTM_his = deval(sol, t_his*0.999 / self.scaleFactor.iv) .* repmat(self.scaleFactor.y, 1, length(t_his));
		end		
        %}
		%% setEventFunc: function description
		function eventFunc = setEventFunc(self, eventPrm)
			self.eventPrm = eventPrm;
			if isfield(eventPrm, 'checkImpact') && eventPrm.checkImpact == 1
				eventFunc = @(t,X) self.detectImpact(t, X, eventPrm.R_imp);
			end
			if isfield(eventPrm, 'checkEscape') && eventPrm.checkEscape == 1
				if ~isfield(eventPrm, 'threshold')
					eventPrm.threshold = 0;
				end
				eventFunc = @(t,X) self.detectEscape(t, X, eventPrm.threshold);
			end
			if isfield(eventPrm, 'checkConvergence') && eventPrm.checkConvergence == 1
				if ~isfield(eventPrm, 'threshold')
					eventPrm.threshold = 1E-05;
					warning(['No threshold value is specified for the convergence check. Instead use ', num2str(eventPrm.threshold)]);
				end
				eventFunc = @(t,X) self.detectConvergence(t, X, eventPrm.f_cmptConvergence, eventPrm.threshold);
			end
		end		
        %{
		%% detectImpact: function description
		function [value, isterminal, direction] = detectImpact(self, t, X, R_imp)
			[rvec, ~] = self.cmpt_rvvec(X, self.mu);
			r = sqrt(sum(rvec.^2)); 
			d = r - R_imp;
			if d < 0
				val_impact = 0;
				if ~isfield(self.eventPrm, 'message') || ~strcmp(self.eventPrm.message, 'suppress')
					warning('impacted')
				end
			else
				val_impact = d;
			end

			value = [val_impact];
			isterminal = [1];
			direction = [0];
		end

		%% detectPlanetApproach: function description
		function [value, isterminal, direction] = detectPlanetApproach(self, t, X, bodyID, r_threshold, terminate)
			t_dim = t * self.scaleFactor.t;
			X_dim = X(1:self.Nx) .* self.scaleFactor.oe;

			[rvec, ~] = self.cmpt_rvvec(X_dim, self.mu);
			rv_planet = self.bodyModels.get_Xpj(t_dim, bodyID);

			value = norm(rvec - rv_planet(1:3), 2) - r_threshold;

			if exist('terminate', 'var')
				isterminal = terminate;
			else
				isterminal = [0];
			end
			direction = [-1];
		end

		%% detectEscape: function description
		function [value, isterminal, direction] = detectEscape(self, t, X, threshold)
			% [rvec, vvec] = self.cmpt_rvvec(X, self.mu);
			% Energy = norm(vvec)^2/2 - self.mu/norm(rvec);
			% if Energy + threshold < 0


			COE = self.cmpt_oe2coe(t, X, self.mu);
			val_escape = 1 - COE(2) - threshold;
			if val_escape<=0
				if ~isfield(self.eventPrm, 'message') || ~strcmp(self.eventPrm.message, 'suppress')
					disp('escaped')
				end
			end

			value = [val_escape];
			isterminal = [1];
			direction = [-1];
		end
        %}
		%% detectConvergence: function description
		function [value, isterminal, direction] = detectConvergence(self, t, X, f_cmptConvergence, threshold)
			t = t * self.scaleFactor.iv;
			X = X.* self.scaleFactor.y;
			maxDiff = f_cmptConvergence(t,X);
			val_converge = maxDiff - threshold;

			if val_converge<=0
				if ~isfield(self.eventPrm, 'message') || ~strcmp(self.eventPrm.message, 'suppress')
					disp('orbit converged')
				end
			end

			value = [val_converge];
			isterminal = [1];
			direction = [-1];
		end
        %{
		%% detectSingularity: function description
		function [value, isterminal, direction] = detectSingularity(self, t, X)
			isterminal = [1];
			direction = [0];
			% check singularity
			if self.singularity.check
				value = self.singularity.fun(t, X);
			else
				value = -1;
			end
			if value >= 0
				disp('dynamics is close to singular')
			end
		end % end for function
        %}
	end
end


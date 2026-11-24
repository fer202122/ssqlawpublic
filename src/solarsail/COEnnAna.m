%% QlawOEnnAna: function description
function output = COEnnAna(t, COE, COEtarg, prm)
% QLAWOENNANA Q-law orbital element double minimization (i.e., nn) analytical
% TAs: true longitude, \in\R^{1\times N_TL}
OrbDiff = COE(1:2) - COEtarg(1:2);
sgm = sign(OrbDiff);
sgm(sgm==0) = 1;

kepvec = [0; 0];

Cprms = prm.Cprms;
alpMinMax = prm.alpMinMax;
g0 = prm.g0;

sma = COE(1);
ecc = COE(2);
inc = COE(3);
omg = COE(5);
Omg = COE(4);
sino = sin(omg);
coso = cos(omg);

tani = tan(inc);
sini = sin(inc);


% Compute Sun vector in Orbital frame (O)
thetas = sunMotion(t, prm.t0jd);
lambda = Omg - thetas;
Rot_ACH2Orb = Rot_body(3,omg)*Rot_body(1,inc)*Rot_body(3,lambda);
% Rot_ACH2Orb = Rot_body(3,omg)*Rot_body(1,inc)*Rot_body(3,Omg);
svec = Rot_ACH2Orb * [1;0;0];

pp = sma*(1-ecc^2);
hh = sqrt(pp*prm.mu);

oenndot = zeros(2,1);
oenndotAn = zeros(2,1);
TA_opt  = zeros(2,1);

% ===analytical f for a===
[oenndot(1), TA_opt(1)] = cmpt_oenndot_a(sgm(1), svec);
oenndotAn(1) = oenndot(1);

%% cmpt_oenndot_a: function description
function [oenndot_a, TA_opt_a] = cmpt_oenndot_a(sgm_a, svec)
	TA_opt_a 	= atan2(sgm_a * svec(1), -sgm_a * svec(2));

	Rot_ACH2LVLH_fst = Rot_body(3,TA_opt_a)*Rot_ACH2Orb;
	r 			= pp / (ecc*cos(TA_opt_a)+1);
	p_oe 		= - sgm_a/hh * 2*sma^2* ([ecc*sin(TA_opt_a), pp/r, 0] * Rot_ACH2LVLH_fst)';
	att_opt 	= attLyapunovFB(p_oe, Cprms);
	% att_opt(1) 	= median([att_opt(1); alpMinMax(:)]);
	acc_opt 	= g0 * work_aSRP_hat(att_opt, Cprms);
	oenndot_a 	= sgm_a*kepvec(1) - dot(p_oe, acc_opt)';
end
% ===analytical f for a===


% ===analytical f for e===
% [edot_nn_num, TA_opt_ecc_num] = cmpt_oenndot_e_num(sgm(2), svec);
% oenndot(2) = edot_nn_num;
% TA_opt(2) = TA_opt_ecc_num;

% [edot_nn, TA_opt_ecc] = cmpt_oenndot_e_num(sgm(2), svec);
[edot_nn, TA_opt_ecc] = cmpt_oenndot_e_ana(sgm(2), svec);
% [edot_nn, TA_opt_ecc] = cmpt_oenndot_e_ana(-1, svec);
% [edot_nn, TA_opt_ecc] = cmpt_oenndot_e_ana(+1, svec); 
oenndot(2) = edot_nn;
TA_opt(2) = TA_opt_ecc;

%% cmpt_oenndot_e_ana: function description
function [oenndot_e, TA_opt_e] = cmpt_oenndot_e_ana(sgm_e, svec)

	% resolve ambiguity by comparing two values
	TA_opt_e0 	= 1/2*atan2(sgm_e*svec(1), -sgm_e*svec(2)) + [0, pi];

	% first true anomaly
	Rot_ACH2LVLH_fst = Rot_body(3,TA_opt_e0(1))*Rot_ACH2Orb;
	r 			= pp / (ecc*cos(TA_opt_e0(1))+1);
	p_oe 		= - sgm_e/hh * ([pp*sin(TA_opt_e0(1)), (pp+r)*cos(TA_opt_e0(1)) + r*ecc, 0] * Rot_ACH2LVLH_fst)';
	att_opt 	= attLyapunovFB(p_oe, Cprms);
	% att_opt(1) 	= median([att_opt(1); alpMinMax(:)]);
	acc_opt 	= g0 * work_aSRP_hat(att_opt, Cprms);
	oenndot_e1 	= sgm_e*kepvec(2) - dot(p_oe, acc_opt)';

	% another true anomaly
	Rot_ACH2LVLH_fst = Rot_body(3,TA_opt_e0(2))*Rot_ACH2Orb;
	r 			= pp / (ecc*cos(TA_opt_e0(2))+1);
	p_oe 		= - sgm_e/hh * ([pp*sin(TA_opt_e0(2)), (pp+r)*cos(TA_opt_e0(2)) + r*ecc, 0] * Rot_ACH2LVLH_fst)';
	att_opt 	= attLyapunovFB(p_oe, Cprms);
	% att_opt(1) 	= median([att_opt(1); alpMinMax(:)]);
	acc_opt 	= g0 * work_aSRP_hat(att_opt, Cprms);
	oenndot_e2 	= sgm_e*kepvec(2) - dot(p_oe, acc_opt)';

	[oenndot_e, idx] = min([oenndot_e1, oenndot_e2]);
	TA_opt_e 	= TA_opt_e0(idx);
end

%% cmpt_oenndot_e_ana_tests: function description
function [oenndot_e, TA_opt_e] = cmpt_oenndot_e_ana_tests(sgm_e, svec)
	TA_opt_e 	= 1/2*atan2(sgm_e*svec(1), -sgm_e*svec(2));

	% BEGIN -- resolve ambiguity based on more positive second-order derivative
	theta = atan2(svec(1), svec(2));
	if sgm_e == 1
		if theta < 0
			n=1;
		else
			n=0;
		end
	else
		n = 0;
	end
	TA_opt_e = TA_opt_e + n * pi;
	% END

	% BEGIN -- resolve ambiguity based on first derivative closer to zero
	% if svec(1)*sin(TA_opt_e + theta) > 0
	% 	TA_opt_e = TA_opt_e + pi;
	% end
	% END


	% %% first-order approx -- with r=p approximation
	% f1 = - sgm_e * norm(svec(1:2),2) * (  sin(2*TA_opt_e + theta) + ecc * sin(TA_opt_e + theta));
	% g1 = - sgm_e * norm(svec(1:2),2) * (2*cos(2*TA_opt_e + theta) + ecc * cos(TA_opt_e + theta));
	% TA_opt_e = TA_opt_e - f1/g1;

	% %% first-order approx -- with e^2=0 approximation
	% f1 = - sgm_e * ((2*sin(2*TA_opt_e + theta) + ecc * (3 + cos(2*TA_opt_e)) * sin(TA_opt_e + theta)));
	% g1 = - sgm_e * ((4*cos(2*TA_opt_e + theta) + ecc * (3 + cos(2*TA_opt_e)) * cos(TA_opt_e + theta) - 2 * ecc * sin(2*TA_opt_e) * sin(TA_opt_e + theta)));
	% f1 = - sgm_e * ((ecc * (3 + cos(2*TA_opt_e)) * sin(TA_opt_e + theta)));
	% g1 = - sgm_e * ((4*cos(2*TA_opt_e + theta) + ecc * (3 + cos(2*TA_opt_e)) * cos(TA_opt_e + theta) - 2 * ecc * sin(2*TA_opt_e) * sin(TA_opt_e + theta)));
	% TA_opt_e = TA_opt_e - f1/g1;

	% %% first-order approx -- with no approximation
	% f1 = - sgm_e * (4*ecc^2*svec(1) + 2 *norm(svec(1:2),2) * (2*sin(2*TA_opt_e + theta) + ecc * (3 + cos(2*TA_opt_e)) * sin(TA_opt_e + theta)));
	% g1 = - sgm_e * (2 *norm(svec(1:2),2) * (4*cos(2*TA_opt_e + theta) + ecc * (3 + cos(2*TA_opt_e)) * cos(TA_opt_e + theta) - 2 * ecc * sin(2*TA_opt_e) * sin(TA_opt_e + theta)));
	% TA_opt_e = TA_opt_e - f1/g1;


	Rot_ACH2LVLH_fst = Rot_body(3,TA_opt_e)*Rot_ACH2Orb;
	r 			= pp / (ecc*cos(TA_opt_e)+1);
	p_oe 		= - sgm_e/hh * ([pp*sin(TA_opt_e), (pp+r)*cos(TA_opt_e) + r*ecc, 0] * Rot_ACH2LVLH_fst)';
	att_opt 	= attLyapunovFB(p_oe, Cprms);
	att_opt(1) 	= median([att_opt(1); alpMinMax(:)]);
	acc_opt 	= g0 * work_aSRP_hat(att_opt, Cprms);
	oenndot_e 	= sgm_e*kepvec(2) - dot(p_oe, acc_opt)';

end

%% cmpt_oenndot_e_num: function description
function [oenndot_e, TA_opt_e] = cmpt_oenndot_e_num(sgm_e, svec)
	% NUMERICAL
	TAs = linspace(0,2*pi,200);
	sTA = sin(TAs);
	cTA = cos(TAs);
	sTh = sin(omg + TAs);
	tinc = tan(inc);
	rs = pp ./ (1+ecc*cTA);


	%% original objective values
	% temp 		= [-rs .*sTA .* (cTA + ecc); pp + rs.*cTA .*(cTA + ecc); zeros(1,numel(TAs)) ];
	% p_oe 		= - sgm_e/hh * Rot_ACH2Orb' * temp;
	% att_opt 	= attLyapunovFB(p_oe, Cprms);
	% att_opt(1,:) = median([att_opt(1,:); repmat(alpMinMax(:), 1, numel(TAs))],1);
	% acc_opt 	= g0 * work_aSRP_hat(att_opt, Cprms);
	% e_temp 		= sgm_e*kepvec(2) - dot(p_oe, acc_opt)';
	%% original objective values

	%% approximate objective values -- based on the dot product between the sunlight and primer vectors
	e_temp = sgm_e*(pp * svec(2) + rs.*(cTA + ecc).*(-svec(1)*sTA + svec(2)*cTA) );
	%% approximate objective values

	[~, idx2opt] = min(e_temp);
	TA_opt_e = TAs(idx2opt);

	Rot_ACH2LVLH_fst = Rot_body(3,TA_opt_e)*Rot_ACH2Orb;
	r 			= pp / (ecc*cos(TA_opt_e)+1);
	p_oe 		= - sgm_e/hh * ([pp*sin(TA_opt_e), (pp+r)*cos(TA_opt_e) + r*ecc, 0] * Rot_ACH2LVLH_fst)';
	att_opt 	= attLyapunovFB(p_oe, Cprms);
	att_opt(1) 	= median([att_opt(1); alpMinMax(:)]);
	acc_opt 	= g0 * work_aSRP_hat(att_opt, Cprms);
	oenndot_e 	= sgm_e*kepvec(2) - dot(p_oe, acc_opt)';
end
% ===analytical f for e===

% output
output.oenndot = oenndot;
output.oenndotAn = oenndotAn;
output.Lopt = TA_opt;

return
% ===analytical f for i===
[oenndot(3), TA_opt(3)] = cmpt_oenndot_i(sgm(3), svec);

%% cmpt_oenndot_i: function description
function [oenndot_i, TA_opt_i] = cmpt_oenndot_i(sgm_i, svec)

	%% OLD ONE
	% if -sgm_i*svec(3) > 0
	% 	TA_opt_i = -omg + asin(-ecc*sino);
	% else
	% 	TA_opt_i = pi - omg - asin(-ecc*sino);
	% end
	%% OLD ONE

	TA_opt_i 	= pi/2 -omg + sign(sgm_i*svec(3)) * (asin(ecc*sino) + pi/2);

	Rot_ACH2LVLH_fst = Rot_body(3,TA_opt_i)*Rot_ACH2Orb;
	r 			= pp / (ecc*cos(TA_opt_i)+1);
	p_oe 		= - sgm_i/hh * ([0, 0, r*cos(TA_opt_i+omg)] * Rot_ACH2LVLH_fst)';
	att_opt 	= attLyapunovFB(p_oe, Cprms);
	att_opt(1) 	= median([att_opt(1); alpMinMax(:)]);
	acc_opt 	= g0 * work_aSRP_hat(att_opt, Cprms);
	oenndot_i 	= sgm_i*kepvec(3) - dot(p_oe, acc_opt)';
end
% ===analytical f for i===


% ===analytical f for Omega===
[oenndot(5), TA_opt(5)] = cmpt_oenndot_O(sgm(5), svec);

%% cmpt_oenndot_O: function description
function [oenndot_O, TA_opt_O] = cmpt_oenndot_O(sgm_O, svec)

	%% OLD ONE
	% if -sgm_O*svec(3)/sini > 0
	% 	TA_opt_O = -omg + acos(-ecc*coso);
	% else
	% 	TA_opt_O = 2*pi - omg - acos(-ecc*coso);
	% end
	%% OLD ONE

	TA_opt_O 	= pi-omg + sign(sgm_O*svec(3)/sini) * acos(ecc*coso);

	Rot_ACH2LVLH_fst = Rot_body(3,TA_opt_O)*Rot_ACH2Orb;
	r 			= pp / (ecc*cos(TA_opt_O)+1);
	p_oe 		= - sgm_O/(hh*sini) * ([0, 0, r*sin(TA_opt_O+omg)] * Rot_ACH2LVLH_fst)';
	att_opt 	= attLyapunovFB(p_oe, Cprms);
	att_opt(1) 	= median([att_opt(1); alpMinMax(:)]);
	acc_opt 	= g0 * work_aSRP_hat(att_opt, Cprms);
	oenndot_O 	= sgm_O*kepvec(5) - dot(p_oe, acc_opt)';
end
% ===analytical f for Omega===



% ===analytical f for omega===
[oenndot(4), TA_opt(4)] = cmpt_oenndot_o_num(sgm(4), svec);

%% cmpt_oenndot_o_ana: function description
function [oenndot_o, TA_opt_o] = cmpt_oenndot_o_ana(sgm_o, svec)

	% in-plane
	TA_opt_o_in = 1/2 * atan2(-sgm_o*svec(2), -sgm_o*svec(1));
	Rot_ACH2LVLH_fst = Rot_body(3,TA_opt_o_in)*Rot_ACH2Orb;
	r 			= pp / (ecc*cos(TA_opt_o_in)+1);
	%p_oe 		= - sgm_o/(hh*ecc) * ([-pp*cos(TA_opt_o_in), (pp+r)*sin(TA_opt_o_in) , - ecc*r*sin(TA_opt_o_in+omg) / tani] * Rot_ACH2LVLH_fst)';
	p_oe 		= - sgm_o/(hh*ecc) * ([-pp*cos(TA_opt_o_in), (pp+r)*sin(TA_opt_o_in) , 0] * Rot_ACH2LVLH_fst)';
	att_opt 	= attLyapunovFB(p_oe, Cprms);
	att_opt(1) 	= median([att_opt(1); alpMinMax(:)]);
	acc_opt 	= g0 * work_aSRP_hat(att_opt, Cprms);
	oenndot_o_in = sgm_o*kepvec(4) - dot(p_oe, acc_opt)';

	% out-plane
	TA_opt_o_out = pi - omg - sign(sgm_o*svec(3)/tani) * acos(ecc*coso);
	Rot_ACH2LVLH_fst = Rot_body(3,TA_opt_o_out)*Rot_ACH2Orb;
	r 			= pp / (ecc*cos(TA_opt_o_out)+1);
	%p_oe 		= - sgm_o/(hh*ecc) * ([-pp*cos(TA_opt_o_out), (pp+r)*sin(TA_opt_o_out) , - ecc*r*sin(TA_opt_o_out+omg) / tani] * Rot_ACH2LVLH_fst)';
	p_oe 		= - sgm_o/(hh) * ([0, 0, - r*sin(TA_opt_o_out+omg) / tani] * Rot_ACH2LVLH_fst)';
	att_opt 	= attLyapunovFB(p_oe, Cprms);
	att_opt(1) 	= median([att_opt(1); alpMinMax(:)]);
	acc_opt 	= g0 * work_aSRP_hat(att_opt, Cprms);
	oenndot_o_out = sgm_o*kepvec(4) - dot(p_oe, acc_opt)';

	% summary
	% oenndot_o 	= (1 - prm.b) * oenndot_o_in + prm.b * oenndot_o_out;
	oenndot_o 	= (oenndot_o_in + prm.b * oenndot_o_out) / (1 + prm.b);
	TA_opt_o 	= nan;

end

%% cmpt_oenndot_o_num: function description
function [oenndot_o, TA_opt_o] = cmpt_oenndot_o_num(sgm_o, svec)

	% NUMERICAL
	TAs = linspace(0,2*pi,200);
	sTA = sin(TAs);
	cTA = cos(TAs);
	sTh = sin(omg + TAs);
	tinc = tan(inc);
	rs = pp ./ (1+ecc*cTA);

	o_temp = sgm_o * rs / hh .* (1/ecc * sTA .* (- svec(1) .* sTA + svec(2) .*cTA ) - sTh ./ tinc * svec(3) );
	% o_temp = sgm_o*(-(pp+rs.*sTA.^2)*svec(1) + rs.*sTA.*cTA * svec(2) - ecc*rs.*sTh/tinc* svec(3) );
	[~, idx4opt] = min(o_temp);
	TA_opt_o = TAs(idx4opt);

	Rot_ACH2LVLH_fst = Rot_body(3,TA_opt_o)*Rot_ACH2Orb;
	r 			= pp / (ecc*cos(TA_opt_o)+1);
	p_oe 		= - sgm_o/(hh*ecc) * ([-pp*cos(TA_opt_o), (pp+r)*sin(TA_opt_o) , - ecc*r*sin(TA_opt_o+omg) / tani] * Rot_ACH2LVLH_fst)';
	att_opt 	= attLyapunovFB(p_oe, Cprms);
	att_opt(1) 	= median([att_opt(1); alpMinMax(:)]);
	acc_opt 	= g0 * work_aSRP_hat(att_opt, Cprms);
	oenndot_o 	= sgm_o*kepvec(4) - dot(p_oe, acc_opt)';


	% plot for checking the number of minima
		% figure(1)
		% hold on
		% p2 = scatter(rad2deg(TA_opt_o), o_temp(idx4opt), 300, 'k^', 'filled');
		% p1 = plot(rad2deg(TAs), o_temp,'-');
		% p1.Color(4) = 0.2;
		% hold off
	% plot for checking the number of minima

end
% ===analytical f for omega===

end % end for QlawOEnnAna




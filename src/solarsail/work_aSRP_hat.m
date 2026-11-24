%% Fsrp: function description
function a_hat = work_aSRP_hat(att, SRP_coefficients)
C1 = SRP_coefficients(1);
C2 = SRP_coefficients(2);
C3 = SRP_coefficients(3);
ca = cos(att(1,:)); sa = sin(att(1,:));
cb = cos(att(2,:)); sb = sin(att(2,:));

a_hat = ca.*[C1.*ca.^2 + C2.*ca + C3; ...
			(C1.*ca + C2).*sa.*sb; ...
			(C1.*ca + C2).*sa.*cb];

end

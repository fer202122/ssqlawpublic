function QLS = rotmatQLS(t, t0jd, COE, coeFrame, obl)
%ROTMATQLS Rotation matrix from Sun (S) to Local orbital frame (L). There
% are two variants, depending on which frame we express the COEs:
% 1) Earth Frame (E)
% 2) Ecliptic Frame (P)
arguments
    t
    t0jd
    COE
    coeFrame = 'Ecliptic';
    obl = 0.0; % If no obliquity is provided set to zero
end

% Health check
if iscolumn(COE)
   COE = COE'; 
end

% Sun's current position in the Ecliptic frame:
thetas = sunMotion(t, t0jd);

% Extract elements
inc  = COE(:,3);
argp = COE(:,5);
raan = COE(:,4);
tra  = COE(:,6);

switch coeFrame
    case 'Ecliptic'
        lambda = raan - thetas;
        % QLS = Rot_body(3,argp+tra)*Rot_body(1,inc)*Rot_body(3,lambda);
        QLS = pagemtimes(Rot_body(3,argp+tra),...
              pagemtimes(Rot_body(1,inc), ...
              Rot_body(3,lambda)));

    case 'Earth'
        % QLS = Rot_body(3,argp+tra)*Rot_body(1,inc)*Rot_body(3,raan)*...
        %       Rot_body(1,obl)*Rot_body(3,-thetas);
        QLS = pagemtimes(Rot_body(3,argp+tra), ...
              pagemtimes(Rot_body(1,inc), ...
              pagemtimes(Rot_body(3,raan), ...
              pagemtimes(Rot_body(1,-obl), ...
                         Rot_body(3,-thetas)))));


    otherwise
        error('Frame not implemented')

end



end


classdef magicformula
    properties
        coeffs
    end

    % Constructor method
    methods (Access = public)
        function obj = magicformula(coeffs)
            obj.coeffs = magicformula.read__(coeffs);
        end
    end

    % Helper methods
    methods (Access = private, Static)
        function c = read__(filename)
            [d,f,ext] = fileparts(filename);
            if strcmpi(ext, ".tir")
                c = simscape.multibody.tirread(filename);
                save([d,f,'.mat'], '-mat', '-struct', 'c')
            elseif strcmpi(ext, ".mat")
                c = load(filename);
                if ~isfield(c, 'LATERAL_COEFFICIENTS')
                    error('Invalid coefficients structure')
                end
            else
                error('Invalid file type. Please provide a .tir or .mat file.');
            end
        end
    end

    % class methods
    methods (Access = public)
        
        function obj = readcoeffs(obj, filename)
            obj.coeffs = obj.read__(filename);
        end

        function [Fy, MzExport] = fy(obj, sa, sr, ia, fz, p)
            % Lateral Force
            % Inputs in SI and radians
            % SAE coordinate system as described on page 62 RCVD
            % sa = slip angle [radians]
            %   positive slip angle gives negative force
            % ia = inclination angle [radians]
            %   positive inclination tilts wheel to the right from behind
            % fz = normal force [Newtons]
            %   Input normal force as negative
            % p = inflation pressure [pascal]
            %   inflation pressure measured in pascal
            
            %% Combine relevant coefficients into one structure
            M = [fieldnames(obj.coeffs.LATERAL_COEFFICIENTS)', ...
                fieldnames(obj.coeffs.SCALING_COEFFICIENTS)', ...
                fieldnames(obj.coeffs.VERTICAL)',...
                fieldnames(obj.coeffs.OPERATING_CONDITIONS);...
                struct2cell(obj.coeffs.LATERAL_COEFFICIENTS)', ...
                struct2cell(obj.coeffs.SCALING_COEFFICIENTS)', ...
                struct2cell(obj.coeffs.VERTICAL)',...
                struct2cell(obj.coeffs.OPERATING_CONDITIONS)'];
            COEFF=struct(M{:});
            
            %% Initializations
            FZ01 = COEFF.LFZO * COEFF.FNOMIN;
            dfz = (abs(fz) - FZ01)/COEFF.FNOMIN*sign(fz);
            dpi = (p - COEFF.NOMPRES)/(COEFF.NOMPRES);
            
            %% Lateral Force Fy
            
            %Pure Slip
            SVyg = fz.*(COEFF.PVY3 + COEFF.PVY4.*dfz).*ia*COEFF.LKYC*COEFF.LMUY;
            SVy0 = fz.*(COEFF.PVY1 + COEFF.PVY2*dfz).*COEFF.LVY*COEFF.LMUY;
            SVy = SVy0 + SVyg;
            Kyg = (COEFF.PKY6 + COEFF.PKY7*dfz)*(1 + COEFF.PPY5*dpi)*fz*COEFF.LKYC;
            Kya = COEFF.PKY1.*COEFF.FNOMIN.*(1 + COEFF.PPY1.*dpi).*sin(COEFF.PKY4.*atan(fz./((COEFF.PKY2 + COEFF.PKY5.*ia.^2)*(1 + COEFF.PPY2.*dpi)*COEFF.FNOMIN))).*(1-COEFF.PKY3.*abs(ia)).*COEFF.LYKA;
            SHy0 = (COEFF.PHY1 + COEFF.PHY2*dfz)*COEFF.LHY;
            SHyg = (Kyg*ia - SVyg)/Kya;
            SHy = SHy0 + SHyg;
            alphay = sa + SHy;

            Ey = (COEFF.PEY1 + COEFF.PEY2*dfz)*(1 + COEFF.PEY5*ia.^2 - (COEFF.PEY3 + COEFF.PEY4.*ia).*sign(alphay)).*COEFF.LEY;
            Mewy = (COEFF.PDY1 + COEFF.PDY2*dfz).*(1 - COEFF.PDY3.*ia.^2).*(1 + COEFF.PPY3.*dpi + COEFF.PPY4.*dpi^2).*COEFF.LMUY;
            Dy = Mewy*fz;
            Cy = COEFF.PCY1*COEFF.LCY;
            By = Kya/(Cy*Dy);
            Fyp = Dy.*sin(Cy*atan(By.*alphay - Ey.*(By.*alphay - atan(By.*alphay)))) + SVy;
            
            % Combined Slip
            SHyk = COEFF.RHY1 + COEFF.RHY2*dfz;
            Eyk = COEFF.REY1 + COEFF.REY2*dfz;
            Cyk = COEFF.RCY1;
            Byk = (COEFF.RBY1 + COEFF.RBY4.*ia.^2).*cos(atan(COEFF.RBY2.*(sa - COEFF.RBY3))).*COEFF.LYKA;
            kappas = sr + SHyk;
            Gyk = cos(Cyk.*atan(Byk.*kappas - Eyk.*(Byk.*kappas - atan(Byk.*kappas))))/cos(Cyk*atan(Byk*SHyk - Eyk*(Byk*SHyk - atan(Byk*SHyk))));
            DVyk = Mewy.*fz.*(COEFF.RVY1 + COEFF.RVY2.*dfz + COEFF.RVY3.*ia).*cos(atan(COEFF.RVY4.*sa));
            SVyk = DVyk.*sin(COEFF.RVY5.*atan(COEFF.RVY6.*sr))*COEFF.LVYKA;
            
            %Disable Combined slip
            % SVyk = 0;
            % Gyk = 1;
            
            % Fy = Lateral Force
            Fy = Gyk*Fyp + SVyk;

            % Export for Mz calculation
            MzExport.SVyg = SVyg;
            MzExport.SVy0 = SVy0;
            MzExport.SVy = SVy;
            MzExport.Kyg = Kyg;
            MzExport.Kya = Kya;
            MzExport.SHy0 = SHy0;
            MzExport.SHyg = SHyg;
            MzExport.SHy = SHy;
            MzExport.By = By;
            MzExport.Cy = Cy;
             
            end

        function [Fx, MzExport] = fx(obj, sa, sr, ia, fz, p)
            % Longitudinal Force
            % MF 6.1 set of equations describing longitudinal force
            % sr = slip ratio [-]
            %   positive slip ratio gives positive force
            % ia = inclination angle [radians]
            %   positive inclination tilts wheel to the right from behind
            % fz = normal force [Newtons]
            %   Input normal force as negative
            % p = inflation pressure [pascal]
            %   inflation pressure measured in pascal

            %% Combine relevant coefficients into one structure
            M = [fieldnames(obj.coeffs.LONGITUDINAL_COEFFICIENTS)', ...
                fieldnames(obj.coeffs.SCALING_COEFFICIENTS)', ...
                fieldnames(obj.coeffs.VERTICAL)',...
                fieldnames(obj.coeffs.OPERATING_CONDITIONS);...
                struct2cell(obj.coeffs.LONGITUDINAL_COEFFICIENTS)', ...
                struct2cell(obj.coeffs.SCALING_COEFFICIENTS)', ...
                struct2cell(obj.coeffs.VERTICAL)',...
                struct2cell(obj.coeffs.OPERATING_CONDITIONS)'];
            COEFF=struct(M{:});
            
            %% Initializations
            FZ01 = COEFF.LFZO * COEFF.FNOMIN;
            dfz = (fz - FZ01)/COEFF.FNOMIN;
            dpi = (p - COEFF.NOMPRES)/(COEFF.NOMPRES);
            
            %% Longitudinal Force Fx
            
            %Pure Slip
            SVx = (COEFF.PVX1 + COEFF.PVX2*dfz)*fz*COEFF.LVX*COEFF.LMX;
            SHx = (COEFF.PHX1 + COEFF.PHX2*dfz)*COEFF.LHX;
            kappax = sr + SHx;
            Kxk = (COEFF.PKX1 + COEFF.PKX2*dfz)*exp(COEFF.PKX3 * dfz)*(1 + COEFF.PPX1*dpi + COEFF.PPX2*dpi^2)*fz*COEFF.LCX;
            Ex  = (COEFF.PEX1 + COEFF.PEX2*dfz + COEFF.PEX3*dfz^2)*(1 - COEFF.PEX4*sign(kappax))*COEFF.LEX;
            Mewx= (COEFF.PDX1 + COEFF.PDX2*dfz)*(1 - COEFF.PDX3*ia^2)*(1 + COEFF.PPX3*dpi + COEFF.PPX4*dpi^2)*COEFF.LMUX;
            Dx  = Mewx*fz;
            Cx  = COEFF.PCX1*COEFF.LCX;
            Bx = Kxk/(Cx*Dx);
            
            % Combined Slip
            SHxa = COEFF.RHX1;
            Exa  = COEFF.REX1 + COEFF.REX2*dfz;
            Cxa  = COEFF.RCX1;
            Bxa  = (COEFF.RBX1 + COEFF.RBX3*ia^2)*cos(atan(COEFF.RBX2*sr))*COEFF.LXAL;
            alphas = sa + SHxa;
            Gxa  = (cos(Cxa*atan(Bxa.*alphas - Exa*(Bxa.*alphas - atan(Bxa.*alphas)))))/(cos(Cxa*atan(Bxa*SHxa - Exa*(Bxa*SHxa - atan(Bxa*SHxa)))));
            
            % Fx = Longitudinal Force
            Fx = (Dx*sin(Cx*atan(Bx.*kappax - Ex.*(Bx.*kappax - atan(Bx.*kappax)))) + SVx)*Gxa;
            
            % Export for Mz Calculation Later
            MzExport.Kxk = Kxk;
            
            end

        function [Mz] = mz(obj, sa, sr, ia, fz, p)
            % Aligning Moment
            % Inputs in SI and radians
            % SAE coordinate system as described on pade 62 RCVD
            
            % sa = slip angle [radians]
            %   positive slip angle gives negative force
            % ia = inclination angle [radians]
            %   positive inclination tilts wheel to the right from behind
            % fz = normal force [Newtons]
            %   Input normal force as negative
            % p = inflation pressure [pascal]
            %   inflation pressure measured in pascal
            
            %% Combine relevant coefficients into one structure
            M = [
                    fieldnames(obj.coeffs.ALIGNING_COEFFICIENTS)', ...
                    fieldnames(obj.coeffs.SCALING_COEFFICIENTS)', ...
                    fieldnames(obj.coeffs.VERTICAL)',...
                    fieldnames(obj.coeffs.OPERATING_CONDITIONS)', ...
                    fieldnames(obj.coeffs.DIMENSION)';...
    
                    struct2cell(obj.coeffs.ALIGNING_COEFFICIENTS)', ...
                    struct2cell(obj.coeffs.SCALING_COEFFICIENTS)', ...
                    struct2cell(obj.coeffs.VERTICAL)',...
                    struct2cell(obj.coeffs.OPERATING_CONDITIONS)', ...
                    struct2cell(obj.coeffs.DIMENSION)'
                ];
            COEFF=struct(M{:});
            
            %% Initializations
            FZ01 = COEFF.LFZO * COEFF.FNOMIN;
            dfz = (fz - FZ01)/COEFF.FNOMIN;
            dpi = (p - COEFF.NOMPRES)/(COEFF.NOMPRES);
            
            %% Lateral Force Fy
            [Fyp0_Gyk0, FYCOEFF] = obj.fy(sa, sr, 0, fz, p);
            
            %% Self-Aligning Moment
            alphaM = sa;
            SHt = COEFF.QHZ1 + COEFF.QHZ2*dfz + (COEFF.QHZ3 + COEFF.QHZ4*dfz)*ia;
            alphat = alphaM + SHt;
            alphar = alphaM + FYCOEFF.SHy + FYCOEFF.SVy/FYCOEFF.Kya;

            %% Pure/Combined Slip
            pure_slip = false;
            if pure_slip
                % Pure Slip
                alphateq = alphat;
                alphareq = alphar;
                s = 0;
                Fx = 0;
            else
                [Fx, FXCOEFF] = obj.fx(sa, sr, ia, fz, p);
                Fy = obj.fy(sa, sr, ia, fz, p);
                % Combined Slip
                alphateq = atan(sqrt((tan(alphat)).^2 + (FXCOEFF.Kxk/FYCOEFF.Kya)^2*sr.^2)).*sin(alphat);
                alphareq = atan(sqrt((tan(alphar)).^2 + (FXCOEFF.Kxk/FYCOEFF.Kya)^2*sr.^2)).*sin(alphar);
                s = (COEFF.SSZ1 + COEFF.SSZ2*(Fy/COEFF.FNOMIN)+(COEFF.SSZ3 + COEFF.SSZ4*dfz)*ia)*COEFF.UNLOADED_RADIUS*COEFF.LS;                
            end
            
            %% Pneumatic Trail t
            Bt = (COEFF.QBZ1 + COEFF.QBZ2*dfz + COEFF.QBZ3*(dfz^2))*(1 + COEFF.QBZ4 + COEFF.QBZ5*abs(ia))*COEFF.LYKA/COEFF.LMY;
            Ct = COEFF.QCZ1;
            Dt = (COEFF.QDZ1 + COEFF.QDZ2*dfz)*(1 - COEFF.PPZ1*dpi)*(1 + COEFF.QDZ3*ia + COEFF.QDZ4*ia^2)*fz*COEFF.UNLOADED_RADIUS/COEFF.FNOMIN*COEFF.LTR;
            Et = (COEFF.QEZ1 + COEFF.QEZ2*dfz + COEFF.QEZ3*dfz^2)*(1 + (COEFF.QEZ4 + COEFF.QEZ5*ia)*(2/pi)*atan(Bt*Ct*alphat));
            t = Dt.*cos(Ct.*atan(Bt.*alphateq - Et.*(Bt.*alphateq - atan(Bt.*alphateq)))).*cos(alphaM);
            
            %% Residual Moment Mzr
            Dr = ((COEFF.QDZ6 + COEFF.QDZ7*dfz)*COEFF.LRES + (COEFF.QDZ8 + COEFF.QDZ9*dfz)*(1 - COEFF.PPZ2*dpi)*ia*COEFF.LKZC + (COEFF.QDZ10 + COEFF.QDZ11*dfz)*ia*abs(ia)*COEFF.LKZC)*fz*COEFF.UNLOADED_RADIUS*COEFF.LMY;
            Br = COEFF.QBZ9*COEFF.LYKA/COEFF.LMY + COEFF.QBZ10*FYCOEFF.By*FYCOEFF.Cy;
            Mzr = Dr.*cos(atan(Br.*alphareq)).*cos(alphaM);
            
            %% Mz
            Mz = -t.*Fyp0_Gyk0 + Mzr + s.*Fx;
            
        end
    end
end
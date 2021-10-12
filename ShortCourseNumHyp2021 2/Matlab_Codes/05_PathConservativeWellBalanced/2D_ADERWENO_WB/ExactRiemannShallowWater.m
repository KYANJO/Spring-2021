% =================================================================== 
% 
%  Exact Riemann solver for the shallow water equations
%  according to 
% 
%  E.F. Toro, Shock-Capturing Methods for Free-Surface Shallow Flows, 
%  Wiley, 2001.
%
%  MATLAB version written by Michael Dumbser 
%
%  please cite the above reference when using this code. 
% 
% =================================================================== 
% 
function [d,u,psi] = ExactRiemannShallowWater(dl,dr,ul,ur,psiL,psiR,g,s)
    %
    % Compute celerity on left and right states
    %
    cl = sqrt(g*dl);
    cr = sqrt(g*dr);
    %
    % Use the "depth positivity condition" to identify
    % type of data and thus of solution and to call
    % appropriate exact solver
    %
    DCRIT = (ur-ul) - 2.0*(cl+cr);
    %
    if( (dl<=0) || (dr<=0) || (DCRIT>=0) )
        disp('Dry bed state generated.')
        return 
    else
        [d,u,psi] = WetBed(s,cl,dl,ul,psiL,cr,dr,ur,psiR,g);
    end
end

function [d,u,psi] = WetBed(s,CL,DL,UL,psiL,CR,DR,UR,psiR,g)
    %
    % Iteration control
    %
    NITER = 1000;
    tol   = 1e-14;
    %
    % Find starting value for iteration
    %
    [DS] = STARTE_SW(CL,DL,UL,CR,DR,UR,g);
    %
    % Store starting value in D0
    %
    D0 = DS;
    %
    % Start iteration
    %
    for IT = [1:NITER]
        [FL,FLD] = GEOFUN(DS,DL,CL,g);
        [FR,FRD] = GEOFUN(DS,DR,CR,g);
        DS  = DS - (FL + FR + UR-UL)/(FLD + FRD);
        CHA = abs(DS-D0)/(0.5*(DS+D0));
        if(CHA<=tol)
            break
        end
        if(DS<0)
            DS = TOL;
        end
        D0 = DS;
    end
    %
    % Converged solution for depth DS in Star Region.
    % Compute velocity US in Star Region
    %
    US = 0.5*(UL + UR) + 0.5*(FR - FL);
    %
    CS = sqrt(g*DS);
    %
    % Sample solution of the wave structure at s=x/t
    %
    [d,u,psi] = SAMWET(s,CL,DL,UL,psiL,CR,DR,UR,psiR,CS,DS,US,g);
    %
end

function [F,FD] = GEOFUN(D,DK,CK,g)
    if(D<=DK)
        % Wave is rarefaction wave (or depression)
        C  = sqrt(g*D);
        F  = 2.0*(C-CK);
        FD = g/C;
    else
        % Wave is shock wave (or bore)
        GES = sqrt(0.5*g*(D+DK)/(D*DK));
        F   = (D-DK)*GES;
        FD  = GES - 0.25*g*(D-DK)/(GES*D*D);
    end
end

function [DS] = STARTE_SW(CL,DL,UL,CR,DR,UR,g)
    DMIN = min(DL,DR);
    % Use Two-Rarefaction (TRRS) solution as starting value
    DS = (1.0/g)*(0.5*(CL+CR)-0.25*(UR-UL))^2;
    if(DS<=DMIN)
        % Use Two-Rarefaction (TSRS) approximation as
        % starting value
    else
        % Use two-shock (TSRS) solution as starting value
        % with DS as computed from TRRS as estimate
        GEL = sqrt(0.5*g*(DS+DL)/(DS*DL));
        GER = sqrt(0.5*g*(DS+DR)/(DS*DR));
        DS  = (GEL*DL + GER*DR - (UR-UL))/(GEL + GER);
    end
end

function [D,U,psi] = SAMWET(S,CL,DL,UL,psiL,CR,DR,UR,psiR,CS,DS,US,g);
    if(S<=US)
        % *****************
        % Sample left wave
        % *****************
        psi=psiL; 
        if(DS>=DL)
            %
            % Left shock
            %
            QL = sqrt((DS + DL)*DS/(2.0*DL*DL));
            SL = UL - CL*QL;
            %
            if(S<=SL)
                % Sample point lies to the left of the shock
                D = DL;
                U = UL;
            else
                % Sample point lies to the right of the shock
                D = DS;
                U = US;
            end
        else
            %
            % Left rarefaction
            %
            SHL = UL - CL;
            %
            if(S<=SHL)
                % Sample point lies to the right of the
                % rarefaction
                D = DL;
                U = UL;
            else
                %
                STL = US - CS;
                %
                if(S<=STL)
                    % Sample point lies inside the rarefaction
                    U = (UL + 2.0*CL + 2.0*S)/3.0;
                    C = (UL + 2.0*CL - S)/3.0;
                    D = C*C/g;
                else
                    % Sample point lies in the STAR region
                    D = DS;
                    U = US;
                end
            end
        end
    else
        % *****************
        % Sample right wave
        % *****************
        psi = psiR; 
        if(DS>=DR)
            %
            % Right shock
            %
            QR = sqrt((DS + DR)*DS/(2.0*DR*DR));
            SR = UR + CR*QR;
            %
            if(S>=SR)
                % Sample point lies to the right of the shock
                D = DR;
                U = UR;
            else
                % Sample point lies to the left of the shock
                D = DS;
                U = US;
            end
        else
            %
            % Right rarefaction
            %
            SHR = UR + CR;
            %
            if(S>=SHR)
                % Sample point lies to the right of the
                % rarefaction
                D = DR;
                U = UR;
            else
                %
                STR = US + CS;
                %
                if(S>=STR)
                    % Sample point lies inside the rarefaction
                    U = (UR  - 2.0*CR + 2.0*S)/3.0;
                    C = (-UR + 2.0*CR + S)/3.0;
                    D = C*C/g;
                else
                    % Sample point lies in the STAR region
                    D = DS;
                    U = US;
                end
            end
        end
    end
end

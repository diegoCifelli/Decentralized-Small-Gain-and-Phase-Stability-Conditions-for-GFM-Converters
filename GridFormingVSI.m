    % This class defines the model of a grid-forming voltage source inverter.

% Author(s): Yitong Li

%% Notes
%
% The model is in load convention.
% 
% The model is in admittance form.
%
% dw means the derivative of w

%% Class

classdef GridFormingVSI < SimplusGT.Class.ModelAdvance
    
    % For temporary use
    properties(Access = protected)
        v_od_r;
        v_oq_r;
        P0;
        Q0;
    end
    
    methods
        % constructor
        function obj = GridFormingVSI(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Static)
        
        function [State,Input,Output] = SignalList(obj)
            State  = {'i_ld','i_lq','x1_pr_d','x2_pr_d','x1_pr_q','x2_pr_q','v_od','v_oq','i_od','i_oq','v_d_ref',...
                        'w','theta', 'p_f', 'q_f', 'i_ld_r', 'i_lq_r'};   
            Input  = {'v_d','v_q','P0'};
            Output = {'i_d','i_q','w','theta'};
        end
        
        % Calculate the equilibrium
        function [x_e,u_e,xi] = Equilibrium(obj)
         	% Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            w   = obj.PowerFlow(5);
            
            % Get parameter
            xwLf    = obj.Para(1);
            Rf      = obj.Para(2);
            xwCf    = obj.Para(3);
            xwLc    = obj.Para(4);
            Rc      = obj.Para(5);
            Xov     = obj.Para(6);
            Rov     = obj.Para(10);
            kp_pr   = 0.5;
            ki_pr   = 30;
            kp_q     = 0.01;
            ki_q     = 1.5;
            xJ      = obj.Para(8);
            xD      = obj.Para(7);
            Q_ctr_ena = obj.Para(9);
            W0      = obj.Para(11);

           
            % Calculate parameters
            Lf = xwLf/W0;
            Cf = xwCf/W0;
            Lc = xwLc/W0;
            
            v_gd = V;
            v_gq = 0;
            i_od = P/V;
            i_oq = -Q/V;
            
            v_gdq = v_gd + 1i*v_gq;
            i_odq = i_od + 1i*i_oq;
            v_odq = v_gdq - i_odq*(Rc + 1i*w*Lc);
            i_cdq = v_odq*(1i*w*Cf);
            i_ldq = i_odq - i_cdq;
            e_dq  = v_odq - i_ldq*(Rf + 1i*w*Lf);
            
            i_ld = -real(i_ldq);
            i_lq = -imag(i_ldq);
            v_od = real(v_odq);
            v_oq = imag(v_odq);
            e_d = real(e_dq);
            e_q = imag(e_dq);
            v_gd = real(v_gdq);
            v_gq = imag(v_gdq);

            theta = xi;            
            obj.P0 = (v_gd*i_od + v_gq*i_oq)*(-1);
            obj.Q0 = (-v_gd*i_oq + v_gq*i_od)*(-1);
            p_f = obj.P0 ;
            q_f = obj.Q0 ;

            % PR current control
            x1 = e_q/ki_pr;
            x2 = e_d/ki_pr;
            error_i_ld = 0;
            error_i_lq = 0;
            x1_pr_d = x1;
            x2_pr_d = x2;
            x1_pr_q = -x2;
            x2_pr_q = x1;

            
            % Virtual admittance
            i_ld_r = i_ld + error_i_ld;
            i_lq_r = i_lq + error_i_lq;
            error_v_od = Rov*i_ld_r - Xov*i_lq_r;
            error_v_oq = Rov*i_lq_r + Xov*i_ld_r;

            % Voltage reference
            v_od_r = v_od + error_v_od;
            v_oq_r = v_oq + error_v_oq;
            obj.v_od_r = v_od_r;
            obj.v_oq_r = v_oq_r;

            % Q PI control
            v_od_r_i = v_od_r/ki_q;
            
            % Notes:
            % Ideally, v_oq_r = 0 should be valid. So, this equilibrium
            % calculation has to be corrected.
            % The P-F and Q-V droop has not been considerred, which should
            % also be corrected.
            
            % Get equilibrium
            x_e = [-i_ld; -i_lq; x1_pr_d; x2_pr_d; x1_pr_q; x2_pr_q; v_od; v_oq; i_od; i_oq; v_od_r_i; w; theta; p_f; q_f; -i_ld_r; -i_lq_r];
            u_e = [v_gd; v_gq; 0];
        end
        
        % State space model
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
        	% Get the power PowerFlow values
            V	= obj.PowerFlow(3);
            
            % Get input
            v_gd   = u(1);
            v_gq   = u(2);
            dPref = u(3);

            % Get state
            i_ld   = x(1);
            i_lq   = x(2);
            x1_pr_d = x(3);
            x2_pr_d = x(4);
            x1_pr_q = x(5);
            x2_pr_q = x(6);
            v_od   = x(7);
            v_oq   = x(8); 
            i_od = x(9);
            i_oq = x(10);
            v_od_r_i = x(11);
            w = x(12);
            theta = x(13);
            p_f = x(14);
            q_f = x(15);
            i_ld_r = x(16);
            i_lq_r = x(17);
            
            % Get parameters
            % Get parameter
            xwLf    = obj.Para(1);
            Rf      = obj.Para(2);
            xwCf    = obj.Para(3);
            xwLc    = obj.Para(4);
            Rc      = obj.Para(5);
            Xov     = obj.Para(6);
            Rov     = obj.Para(10);
            kp_pr   = 0.5;
            ki_pr   = 30;
            kp_q     = 0.01;
            ki_q     = 5;
            xJ      = obj.Para(8);
            xD      = obj.Para(7);
            Q_ctr_ena = obj.Para(9);
            W0      = obj.Para(11);

            
            % Update paramters
            Lf = xwLf/W0;
            Cf = xwCf/W0;
            Lc = xwLc/W0;
            J = xJ/W0;   
            D = xD; 
            
            % Saturation setting
            EnableSaturation = 0;
            
            % Frequency limit and saturation
            w_limit_H = W0*1.1;
            w_limit_L = W0*0.9;
            % Capacitor voltage limit
            v_od_limit_H = 1.5;
            v_od_limit_L = -1.5;
            v_oq_limit_H = 1.5;
            v_oq_limit_L = -1.5;    
            % Current reference limit
            i_ld_limit = 1.5;
            i_lq_limit = 1.5;
            % Ac voltage limit
            e_d_limit_H = 1.5;
            e_d_limit_L = -1.5;
            e_q_limit_H = 1.5;
            e_q_limit_L = -1.5;
       
  
            % State space equations
         	% dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1    
            % ### Call state equation: dx/dt = f(x,u)
            
                % Power measurement
                p = (v_gd*i_od + v_gq*i_oq)*(-1);   % (-1) appears because the model is in load convention
                q = (-v_gd*i_oq + v_gq*i_od)*(-1);

                % Virtual Synchonoues Machine /VSM)
                P0     = obj.P0+dPref;
                Q0     = obj.Q0;
                dp_f = 0;
                dw = ((W0-w)*D + (P0 - p))/2/J;   
                
                % Limitation for w
                if EnableSaturation
                     if (w >= w_limit_H && dw >=0 ) || (w <= w_limit_L && dw <= 0 )
                        dw = 0;                  
                     end
                end
                
                % Reactive power / Voltage Controller
                switch Q_ctr_ena      
                    case 1
                        dq_f = 0;
                        error_Q = Q0-q;
                        v_od_r = kp_q*error_Q + ki_q*v_od_r_i;                  % Q-V droop with LPF
                        dv_od_r_i = error_Q;
                        % Limitation for v_od_r
                        if EnableSaturation
                            v_od_r = min(v_od_r,v_od_limit_H);
                            v_od_r = max(v_od_r,v_od_limit_L);
                        end  
                    case 0
                        dq_f = 0;
                        v_od_r = obj.v_od_r;                                % No Q-V droop
                        dv_od_r_i = 0;
                        % Limitation for v_od_r
                        if EnableSaturation
                            v_od_r = min(v_od_r,v_od_limit_H);
                            v_od_r = max(v_od_r,v_od_limit_L);
                        end                        
                    otherwise
                        error(['Error'])
                end

                % Limitation for v_oq_r
                if EnableSaturation
                    v_oq_r = min(v_oq_r,v_oq_limit_H);
                    v_oq_r = max(v_oq_r,v_oq_limit_L);
                end  
                                
                % Virtual admittance
                error_v_od = v_od_r - v_od;
          	    error_v_oq = obj.v_oq_r - v_oq; 
                wctrl = w;
                di_ld_r = -(error_v_od - Rov*-i_ld_r + wctrl*Xov/W0*-i_lq_r)/(Xov/W0);
                di_lq_r = -(error_v_oq - Rov*-i_lq_r - wctrl*Xov/W0*-i_ld_r)/(Xov/W0);

                % Current saturation
                if EnableSaturation
                    i_ld_r = min(i_ld_r,i_ld_limit);
                    i_ld_r = max(i_ld_r,-i_ld_limit);
                    i_lq_r = min(i_lq_r,i_lq_limit);
                    i_lq_r = max(i_lq_r,-i_lq_limit);
                end

                % AC current control (PR)
                wo_pr = W0;
                zeta_pr = 0;
                error_i_ld = i_ld_r-i_ld;
                error_i_lq = i_lq_r-i_lq;
                dx1_pr_d = wo_pr*x2_pr_d + wctrl*x1_pr_q;
                dx2_pr_d = -wo_pr*x1_pr_d - 2*zeta_pr*wo_pr*x2_pr_d  + wctrl*x2_pr_q + error_i_ld*(-1);
                dx1_pr_q = -wctrl*x1_pr_d + wo_pr*x2_pr_q;
                dx2_pr_q = -wctrl*x2_pr_d - wo_pr*x1_pr_q - 2*zeta_pr*wo_pr*x2_pr_q + error_i_lq*(-1);
                e_d = ki_pr*x2_pr_d + kp_pr*error_i_ld*(-1);
                e_q = ki_pr*x2_pr_q + kp_pr*error_i_lq*(-1);


                % Ac voltage (duty cycle) saturation
                if EnableSaturation
                    e_d = min(e_d,e_d_limit_H);
                    e_d = max(e_d,e_d_limit_L);
                    e_q = min(e_q,e_q_limit_H);
                    e_q = max(e_q,e_q_limit_L);
                end

                % Lf equation
                % e_d - v_od = -(di_ld/dt*Lf + Rf*i_ld - w*Lf*i_lq)
                % e_q - v_oq = -(di_lq/dt*Lf + Rf*i_lq + w*Lf*i_ld)
                di_ld = (v_od - e_d - Rf*i_ld + wctrl*Lf*i_lq)/Lf;
                di_lq = (v_oq - e_q - Rf*i_lq - wctrl*Lf*i_ld)/Lf;

                % Cf equation - MOVED EXTERNALALY
                % -(i_ld - i_od) = Cf*dv_cd/dt - w*Cf*v_cq
                % -(i_lq - i_oq) = Cf*dv_cq/dt + w*Cf*v_cd
                dv_od = (-(i_ld - i_od) + wctrl*Cf*v_oq )/Cf;
                dv_oq = (-(i_lq - i_oq) - wctrl*Cf*v_od )/Cf;

                % Lc equation - MOVED EXTERNALALY
                % v_od - v_d = -(Lc*di_od/dt + Rc*i_od - w*Lc*i_oq)
                % v_oq - v_q = -(Lc*di_oq/dt + Rc*i_oq + w*Lc*i_od)
                di_od = (v_gd - v_od - Rc*i_od + wctrl*Lc*i_oq)/Lc;
                di_oq = (v_gq - v_oq - Rc*i_oq - wctrl*Lc*i_od)/Lc;


                % Phase angle
                dtheta = w;  

                % dx              
                f_xu = [di_ld; di_lq; dx1_pr_d; dx2_pr_d;dx1_pr_q; dx2_pr_q; dv_od; dv_oq; di_od; di_oq; dv_od_r_i; dw; dtheta; dp_f; dq_f; di_ld_r; di_lq_r];
                Output = f_xu;
                
            elseif CallFlag == 2     
            % ### Call output equations: y = g(x,u)
                g_xu = [i_od; i_oq; w; theta];
                Output = g_xu;
            end
            
        end
        
    end
end
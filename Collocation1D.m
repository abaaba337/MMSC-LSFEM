


classdef Collocation1D < handle 

    properties                     
        dep_num  ;      % m, number of dependent variables
        num_eq   ;      % number of equations
        num_cond ;      % number of conditions
        A1 ; A0 ;  % coefficient matrices
        B1 ; B0 ;  % coefficient matrices for boundary
        f  ; g  ;       % force term and boundary term
        Omega ;         % Omega = [a,b]
       
        u ;             % true solution and its derivatives (if given). 
        Du ;

        basis_num ;   % n, number of global basis
        ele_num   ;   % number of elements
        gm ;
        C  ;          % Solution coefficient matrix of size dep_num by basis_num

        Omega_assemble ; % Omega_assemble.B : B_Omega, Omega_assemble.d : d_Omega, Omega_assemble.N : N_Omega
        Gamma_assemble ; % Gamma_assemble.B : B_Gamma, Gamma_assemble.d : d_Gamma, Gamma_assemble.N : N_Gamma

    end



    methods

        function obj = Collocation1D(dep_num, Omega, ...
                                     A1, A0, f, num_eq, ...
                                     B1, B0, g, num_cond, ...
                                     u, Du) 
            
            arguments 
                dep_num
                Omega

                A1 = []
                A0 = []
                f = [] 
                num_eq = []
                B1 = []
                B0 = []
                g = []
                num_cond = []
                
                u = []
                Du = []
            end
            
            obj.dep_num = dep_num ;
            obj.A1 = A1 ; obj.A0 = A0 ; 
            obj.B1 = B1 ; obj.B0 = B0 ; 
            obj.f = f   ; obj.g = g ;

            if isempty(num_eq)
                obj.num_eq = dep_num ;
            else
                obj.num_eq = num_eq ; 
            end

            if isempty(num_cond)
                obj.num_cond = dep_num ;
            else
                obj.num_cond = num_cond ;
            end

            obj.u = u ;
            obj.Du = Du ;

            obj.Omega = Omega ;
        end

        % This function fit interior conditions after model construction
        function fitInterior(obj, A1, A0, f, num_eq)
           arguments 
                obj
                A1 
                A0 
                f 
                num_eq = obj.dep_num
            end
            obj.A1 = A1 ; obj.A0 = A0 ; obj.f = f;
            obj.num_eq = num_eq ;
        end

        % This function fit boundary conditions after model construction
        function fitBoundary(obj, B1, B0, g, num_cond)
            arguments 
                obj
                B1 
                B0 
                g 
                num_cond = obj.dep_num
            end

            obj.num_cond = num_cond;

            if isempty(B1)
                obj.B1 = @(x,y) zeros(obj.num_cond, obj.dep_num);
            else
                obj.B1 = B1; 
            end

            if isempty(B0)
                obj.B0 = @(x,y) eye(obj.num_cond);
            else
                obj.B0 = B0; 
            end
            
            obj.g = g ; 
        end

        % This function fit solution benchmark after model construction
        function fitTrueSol(obj, u, Du)
            arguments 
                obj
                u
                Du = []
            end

            obj.u = u ; obj.Du = Du ;
        end

        
        % This function construct a simple uniform mesh for the domian CG1.
        function linear_discretize(obj, n)

            obj.gm.Mesh.Nodes = linspace(obj.Omega(1), obj.Omega(2), n) ;
            obj.gm.Mesh.Elements = [ 1:n-1 ; 2:n ] ;
            obj.gm.Mesh.MaxElementSize = abs(obj.Omega(1) - obj.Omega(2)) / (n-1) ;
            obj.gm.Mesh.MinElementSize = abs(obj.Omega(1) - obj.Omega(2)) / (n-1) ;
            obj.gm.Mesh.GeometricOrder = 'linear' ;
            obj.basis_num = n ; 
            obj.ele_num = n-1 ;
        
        end

 
        % This function evaluates B_Omega locally.
        function res = B_loc(obj, x, eleID, position, Lambda)
            % if position == 0 evaluate interior; if position == 1 evaluate
            % boundary
            
            arguments 
                obj
                x              
                eleID              
                position = 0       
                Lambda = eye(obj.num_cond)
            end

            v_globID = obj.gm.Mesh.Elements(:, eleID);
            v = obj.gm.Mesh.Nodes(v_globID) ;
            l_loc  = v(end) - v(1) ;
       
            if position == 0
                if strcmpi(obj.gm.Mesh.GeometricOrder, 'linear')
                    B1_loc = (- obj.A1(x) + (v(2) - x) * obj.A0(x) )/l_loc ;
                    B2_loc = (  obj.A1(x) + (x - v(1)) * obj.A0(x) )/l_loc ;
                end
            elseif position == 1
                if strcmpi(obj.gm.Mesh.GeometricOrder, 'linear')
                    B1_loc = Lambda * (- obj.B1(x) + (v(2) - x) * obj.B0(x) )/l_loc ;
                    B2_loc = Lambda * (  obj.B1(x) + (x - v(1)) * obj.B0(x) )/l_loc ;
                end
            end 
            res = [ B1_loc , B2_loc ] ;
        end

 
        % This function evaluates and assembles global B_Omega.
        function [B_Omega, d_Omega, N_Omega_modified] = Assemble_B_Omega(obj, interior_points, save_result)

           arguments 
               obj
               interior_points
               save_result = true
           end

            n = obj.basis_num ;           % number of global basis
            m = obj.dep_num ;             % number of dependent variables
            N_Omega = length(interior_points) ;  % number of interior points
            mesh = obj.gm.Mesh ; 
        
            Cols = cell(1, N_Omega)    ; Rows = cell(1, N_Omega) ;
            B_Omega = cell(1, N_Omega) ; d_Omega = cell(1, N_Omega) ;
        
            % Loop over all interior points and calculate local B_Omega, d_Omega in
            % parallel
        
            hbar = parfor_progressbar_v1(N_Omega,'Local B_{\Omega} calculation, please wait...') ; % progress bar
            parfor i = 1 : N_Omega
                x = interior_points(i) ;
                eleID = obj.findElement(x, mesh) ;
                iota = mesh.Elements(:, eleID) ; 
                B_Omega{i} = obj.B_loc(x, eleID, 0) ;
                d_Omega{i} = obj.f(x) ;
        
                row_ID = ones(1,obj.num_eq) ;
                col_ID = arrayfun(@(i)(i-1)*m+1 : i*m, iota, 'UniformOutput', false) ;
                col_ID = cat(2, col_ID{:}) ;
                [Cols{i}, Rows{i}] = meshgrid(col_ID, row_ID) ;

                hbar.iterate(1) ;
            end
            close(hbar); % close the progress bar
        
            % Assemble B_Omega, d_Omega
            Rows = cumsum(cat(1, Rows{:})) ; maxRow = max(Rows(:,1)) ; 
            Cols = cat(1, Cols{:}) ;
            B_Omega = cat(1, B_Omega{:}) ; d_Omega = cat(1, d_Omega{:}) ;
        
            Rows = reshape(Rows, [], 1) ; Cols = reshape(Cols, [], 1) ;
            B_Omega = reshape(B_Omega, [], 1) ;
            B_Omega = sparse(Rows, Cols, B_Omega, maxRow, m*n) ;
            N_Omega_modified = round(maxRow/m) ;

            if save_result
                obj.Omega_assemble.B = B_Omega;
                obj.Omega_assemble.d = d_Omega;
                obj.Omega_assemble.N = N_Omega_modified;
            end
        
        end


        % This function evaluates and assembles global B_Gamma.
        function [B_Gamma, d_Gamma, N_Gamma_modified] = Assemble_B_Gamma(obj, boundary_points, Lambda, save_result)
        
            % Lambda m ny m square(weights) for BC

           arguments 
               obj
               boundary_points
               Lambda
               save_result = true
           end
        
            n = obj.basis_num ;           % number of global basis
            m = obj.dep_num ;             % number of dependent variables
            N_Gamma = length(boundary_points) ;   % number of boundary points
            mesh = obj.gm.Mesh ; 
        
            Cols = cell(1, N_Gamma)    ; Rows = cell(1, N_Gamma) ;
            B_Gamma = cell(1, N_Gamma) ; d_Gamma = cell(1, N_Gamma) ;
            
            % Loop over all boundary points and calculate local B_Gamma, d_Gamma in
            % parallel
        
            hbar = parfor_progressbar_v1(N_Gamma,'Local B_{\Gamma} calculation, please wait...') ; % progress bar
            parfor i = 1 : N_Gamma
                x = boundary_points(i) ;
                eleID = obj.findElement(x, mesh) ;
        
                iota = mesh.Elements(:, eleID) ; 
                B_Gamma{i} = obj.B_loc(x, eleID, 1, Lambda) ;
                d_Gamma{i} = Lambda*obj.g(x) ;
        
                row_ID = ones(1,obj.num_cond) ;
                col_ID = arrayfun(@(i)(i-1)*m+1 : i*m, iota, 'UniformOutput', false) ;
                col_ID = cat(2, col_ID{:}) ;
                [Cols{i}, Rows{i}] = meshgrid(col_ID, row_ID) ;
        
                hbar.iterate(1) ;
            end
            close(hbar); % close the progress bar
        
            % Assemble B_Omega, d_Omega
            Rows = cumsum(cat(1, Rows{:})) ; maxRow = max(Rows(:,1)) ; 
            Cols = cat(1, Cols{:}) ;
            B_Gamma = cat(1, B_Gamma{:}) ; d_Gamma = cat(1, d_Gamma{:}) ;
        
            Cols = reshape(Cols, [], 1) ;
            Rows = reshape(Rows, [], 1) ;
            B_Gamma = reshape(B_Gamma, [], 1) ;
            B_Gamma = sparse(Rows, Cols, B_Gamma, maxRow, m*n) ;
            N_Gamma_modified = round(maxRow/m) ;

            if save_result
                obj.Gamma_assemble.B = B_Gamma;
                obj.Gamma_assemble.d = d_Gamma;
                obj.Gamma_assemble.N = N_Gamma_modified;
            end
        
        end


        % This function evaluate the point-wise numerical value according to the solution coefficient matrix C.
        function u_fit = evaluate_LSFEM_elewise(obj, x, C, derivative_order)
            % u_fit is a column vector.

            arguments 
                obj
                x
                C = obj.C
                derivative_order = 0
            end

            if isempty(C)
                C = obj.C;
            end

            mesh = obj.gm.Mesh ;
            eleID = obj.findElement(x, mesh) ;
            v_globID = mesh.Elements(:,eleID) ;
            v = mesh.Nodes(v_globID) ;
            l_loc  = v(end) - v(1) ;

            if strcmpi(obj.gm.Mesh.GeometricOrder, 'linear')
                if derivative_order == 0
                    u_fit = C(:,v_globID(1)) * (v(2) - x)/l_loc + C(:,v_globID(2)) * (x - v(1))/l_loc ;
                elseif derivative_order == 1
                    u_fit = - C(:,v_globID(1)) ./ l_loc + C(:,v_globID(2)) ./ l_loc ;
                end
            end
        
        end


        % This function evaluate the vectorized numerical value according to the solution coefficient matrix C.
        function u_fit = evaluate_LSFEM(obj, X, C, derivative_order, showbar)
            % X vector
            arguments
                obj
                X
                C = obj.C
                derivative_order = 0
                showbar = true
            end

            if isempty(C)
                C = obj.C;
            end

            u_fit = cell(1, length(X)) ;

            if showbar
                hbar = parfor_progressbar_v1(length(X), 'Evaluation, please wait...') ;
                parfor i = 1 : length(X)
                    u_fit{i} = obj.evaluate_LSFEM_elewise(X(i), C, derivative_order) ;
                    hbar.iterate(1) ; 
                end
                close(hbar) ;
            else
                parfor i = 1 : length(X)
                    u_fit{i} = obj.evaluate_LSFEM_elewise(X(i), C, derivative_order) ;
                end
            end
           u_fit = cat(2, u_fit{:}) ;
        end


        % This function evaluate the absolute error || u_true(x,y) - u_hat(x,y) ||2 at a series of fixed points.
        function err = abs_error(obj, x, C, derivative_order)
            arguments
                obj
                x
                C = obj.C
                derivative_order = 0
            end

            if isempty(C)
                C = obj.C;
            end

            [sz1,sz2] = size(x);
            x = reshape(x, [], 1);

            if derivative_order == 0
                err = obj.evaluate_LSFEM(x, C, derivative_order, false) - obj.u(x'); % obj.u() input is a row vector
            elseif derivative_order == 1
                err = obj.evaluate_LSFEM(x, C, derivative_order, false) - obj.Du(x'); % obj.u() input is a row vector
            end
            
            err = sqrt(sum(err .* err, 1))';
            err(isnan(err)) = 0;
            err = reshape(err, sz1, sz2);
        end
    
        
        % This function approximate the L2 error || u_true - u_hat ||L2.
        function err = Hk_error(obj, C, type)
            arguments
                obj
                C = obj.C
                type = 0 % Hk error, type = 0 or 1 for k = 0 or 1
            end

           if isempty(C)
                C = obj.C;
            end

            if type == 0
                err = sqrt(integral(@(x) obj.abs_error(x,C).^2, obj.Omega(1), obj.Omega(2), 'AbsTol', 1e-12));

            elseif type == 1
                err = sqrt(integral(@(x) obj.abs_error(x,C).^2 + obj.abs_error(x,C,1).^2, obj.Omega(1), obj.Omega(2), 'AbsTol', 1e-12));

            end
        end


         % This function reshapes the vectorized c to matrix version C.
         function C = reshape2mat(obj, c)
            C = reshape(c, obj.dep_num, obj.basis_num) ;
         end

    end



    methods (Static)
       
         % This function finds the element (ID) where the point xy is located.
         function eleID = findElement(x, mesh)
             v = mesh.Nodes ;
             n = length(v) ;
             eleID = min(sum(v <= x), n-1) ; 
         end


     end


end
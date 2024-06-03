


classdef Collocation1D < handle 

    properties                     
        dep_num  ;      % m, number of dependent variables
        A1 ; A2 ; A0 ;  % coefficient matrices
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

        function obj = Collocation1D(dep_num, A1, A0, f, g, Omega, u, Du) 
            
            arguments 
                dep_num
                A1
                A0
                f
                g
                Omega
                u = []
                Du = []
            end
            
            obj.dep_num = dep_num ;
            obj.A1 = A1 ; obj.A0 = A0 ; 
            obj.f = f ; obj.g = g ;

            obj.u = u ;
            obj.Du = Du ;

            obj.Omega = Omega ;
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
        function res = B_Omega_loc(obj, x, eleID)

            v_globID = obj.gm.Mesh.Elements(:, eleID);
            v = obj.gm.Mesh.Nodes(v_globID) ;
            l_loc  = v(end) - v(1) ;
       
            if strcmpi(obj.gm.Mesh.GeometricOrder, 'linear')
                B1_loc = (- obj.A1(x) + (v(2) - x) * obj.A0(x) )/l_loc ;
                B2_loc = (  obj.A1(x) + (x - v(1)) * obj.A0(x) )/l_loc ;
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
                B_Omega{i} = obj.B_Omega_loc(x, eleID) ;
                d_Omega{i} = obj.f(x) ;
        
                row_ID = ones(1,m) ;
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


        % This function evaluates B_Gamma locally.
        function res = B_Gamma_loc(obj, x, eleID, Lambda)
            
            v_globID = obj.gm.Mesh.Elements(:,eleID) ;
            v = obj.gm.Mesh.Nodes(v_globID) ;
            
            l_loc  = v(end) - v(1) ;
            
            if strcmpi(obj.gm.Mesh.GeometricOrder, 'linear')
                B1_loc = (v(2) - x)/l_loc * Lambda ;
                B2_loc = (x - v(1))/l_loc * Lambda ;
            end

            res = [ B1_loc , B2_loc ] ;
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
                B_Gamma{i} = obj.B_Gamma_loc(x, eleID, Lambda) ;
                d_Gamma{i} = Lambda*obj.g(x) ;
        
                row_ID = ones(1,m) ;
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
        function u_fit = evaluate_LSFEM_elewise(obj, x, C)
            % u_fit is a column vector.

            arguments 
                obj
                x
                C = obj.C
            end

            mesh = obj.gm.Mesh ;
            eleID = obj.findElement(x, mesh) ;
            v_globID = mesh.Elements(:,eleID) ;
            v = mesh.Nodes(v_globID) ;
            l_loc  = v(end) - v(1) ;

            if strcmpi(obj.gm.Mesh.GeometricOrder, 'linear')
                u_fit = C(:,v_globID(1)) * (v(2) - x)/l_loc + C(:,v_globID(2)) * (x - v(1))/l_loc ;
            end
        
        end


        % This function evaluate the vectorized numerical value according to the solution coefficient matrix C.
        function u_fit = evaluate_LSFEM(obj, X, C)
            % XY = [ X , Y ]  
            arguments
                obj
                X
                C = obj.C
            end

            u_fit = cell(1, length(X)) ;

            hbar = parfor_progressbar_v1(length(X), 'Evaluation, please wait...') ;
            parfor i = 1 : length(X)
                u_fit{i} = obj.evaluate_LSFEM_elewise(X(i), C) ;
                hbar.iterate(1) ; 
            end
            close(hbar) ;

           u_fit = cat(2, u_fit{:}) ;
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



classdef Collocation2D < handle 

    properties                     
        dep_num  ;      % m, number of dependent variables
        num_eq   ;      % number of equations
        num_cond ;      % number of conditions
        A1 ; A2 ; A0 ;  % coefficient matrices for interior
        B1 ; B2 ; B0 ;  % coefficient matrices for boundary
        f  ; g  ;       % force term and boundary term
        Omega ;         % polygonal domain Omega = polyshape({Gamma1X, ...}, {Gamma1Y, ...})
        BoundaryNodeIDs ;
        Edges ;
       
        u ;             % true solution and its derivatives (if given). 
        Dx_u ;
        Dy_u ;

        basis_num ; % n, number of global basis
        ele_num   ;   % number of elements
        gm ;
        C  ;        % Solution coefficient matrix of size dep_num by basis_num

        Omega_assemble ; % Omega_assemble.B : B_Omega, Omega_assemble.d : d_Omega, Omega_assemble.N : N_Omega
        Gamma_assemble ; % Gamma_assemble.B : B_Gamma, Gamma_assemble.d : d_Gamma, Gamma_assemble.N : N_Gamma

    end



    methods

        function obj = Collocation2D(dep_num, Omega, ...
                                     A1, A2, A0, f, num_eq,...
                                     B1, B2, B0, g, num_cond,...
                                     u, Dx_u, Dy_u) 
            
            arguments 
                dep_num
                Omega
                A1 = []
                A2 = []
                A0 = []
                f = [] 
                num_eq = []
                B1 = []
                B2 = []
                B0 = []
                g = []
                num_cond = []
                u = []
                Dx_u = []
                Dy_u = []
            end
            
            obj.dep_num = dep_num ;
            obj.A1 = A1 ; obj.A2 = A2 ; obj.A0 = A0 ; 
            obj.B1 = B1 ; obj.B2 = B2 ; obj.B0 = B0 ; 
            obj.f = f ; obj.g = g ;

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
            obj.Dx_u = Dx_u ;
            obj.Dy_u = Dy_u ;
       
            obj.Omega = Omega ;
            tr = triangulation(obj.Omega) ; 
            obj.gm = fegeometry(tr) ; 
        end

        % This function fit interior conditions after model construction
        function fitInterior(obj, A1, A2, A0, f, num_eq)
           arguments 
                obj
                A1 
                A2 
                A0 
                f 
                num_eq = obj.dep_num
            end
            obj.A1 = A1 ; obj.A2 = A2 ; obj.A0 = A0 ; obj.f = f;
            obj.num_eq = num_eq ;
        end

        % This function fit boundary conditions after model construction
        function fitBoundary(obj, B1, B2, B0, g, num_cond)
            arguments 
                obj
                B1 
                B2 
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

            if isempty(B2)
                obj.B2 = @(x,y) zeros(obj.num_cond, obj.dep_num);
            else
                obj.B2 = B2; 
            end

            if isempty(B0)
                obj.B0 = @(x,y) eye(obj.num_cond);
            else
                obj.B0 = B0; 
            end
            
            obj.g = g ; 
        end

        % This function fit solution benchmark after model construction
        function fitTrueSol(obj, u, Dx_u, Dy_u)
            arguments 
                obj
                u
                Dx_u = []
                Dy_u = []
            end

            
            obj.u = u ; obj.Dx_u = Dx_u ; obj.Dy_u = Dy_u ;

        end


        % This function construct a simple uniform triangular mesh for the domian.
        function discretize(obj, hmax, beta, PolyOrder, update_Omega)

            arguments 
                obj
                hmax                  % target maximum mesh edge length
                beta                  % hmin=hmax/beta target minimum mesh edge length
                PolyOrder = 'linear'  % or 'quadratic'
                update_Omega = true
            end

            hmin=hmax/beta ;
            obj.gm = generateMesh(obj.gm, GeometricOrder=PolyOrder, Hmax=hmax, Hmin=hmin) ; % linear discretization
            obj.basis_num = size(obj.gm.Mesh.Nodes, 2)  ; 
            obj.ele_num = size(obj.gm.Mesh.Elements, 2) ;
            Elements = obj.gm.Mesh.Elements(:);
            counts = groupcounts(Elements);

            % find boundary nodes
            obj.BoundaryNodeIDs = find(counts<=4)';
            id_diff = [(obj.BoundaryNodeIDs(2:end) - obj.BoundaryNodeIDs(1:end-1)), 0];
            end_id = sum(id_diff == 1) + 1;
            obj.BoundaryNodeIDs = obj.BoundaryNodeIDs(1:end_id);
            obj.Edges = obj.boundary_decompose();
            if update_Omega
                Ex = cell2mat(cellfun(@(e) [e(1:end,1); nan]', obj.Edges, 'UniformOutput', false));
                Ey = cell2mat(cellfun(@(e) [e(1:end,2); nan]', obj.Edges, 'UniformOutput', false));
                obj.Omega = polyshape(Ex, Ey); % update Omega
            end
        end


        % This function seperate the boundary nodes from the given mesh, only applicable when there is no holes in the domain.
        function [Mesh_Boundary_Nodes_ID, Mesh_interior_Nodes_ID] = seperate_mesh_nodes(obj, s)

            arguments 
                obj   
                s = 0  % shrink factor. If s = 0, returns the convexhull.
            end
        
            XY = obj.gm.Mesh.Nodes' ;
            Mesh_Boundary_Nodes_ID = boundary(XY,s) ;     % find the bounday nodes of the mesh ID(end) = ID(1)
            is_interior = ~ismember(1:obj.basis_num, Mesh_Boundary_Nodes_ID) ;
            Mesh_interior_Nodes_ID = find(is_interior)' ; % find the interior nodes of the mesh
        
        end


        % This function decompose the boundary into seperate edges.
        function edge = boundary_decompose(obj)
            XY = obj.gm.Mesh.Nodes(:, obj.BoundaryNodeIDs)';
            edge_ID = nearestEdge(obj.gm, XY);
            edge = cell(1, obj.gm.NumEdges) ;
            for i = 1:obj.gm.NumEdges
                edge_node_ID = obj.BoundaryNodeIDs(edge_ID == i);
                edge_node_coord = obj.gm.Mesh.Nodes(:,edge_node_ID)';

                if norm(edge_node_coord(end,:)-edge_node_coord(1,:)) <= obj.gm.Mesh.MaxElementSize
                    edge{i} = [edge_node_coord;edge_node_coord(1,:)];
                else
                    edge{i} = edge_node_coord;
                end
            end
        end

    
        % This function calculate the outer normal vector as a row vector for a given boundary node.
        function n = boundary_normal(obj, xy, position)
            arguments 
                obj
                xy
                position = [1,2]
            end
            edgeID = nearestEdge(obj.gm, xy);
            edge = obj.Edges{edgeID};
            vec = edge(2:end,:) - edge(1:end-1,:);
            ks = vec(:,2)./vec(:,1);
            bs = edge(1:end-1,2) - ks .* edge(1:end-1,1);
            segID = abs(ks * xy(1) + bs - xy(2)) < 1e-6;
            vec = sum(vec(segID,:), 1);
            n = [vec(2), -vec(1)]/norm(vec);
            n = n(position);
        end


        % This function evaluates B_Omega and B_Gamma locally.
        function res = B_loc(obj, xy, eleID, position, Lambda)
            % if position == 0 evaluate interior; if position == 1 evaluate
            % boundary
            
            arguments 
                obj
                xy                
                eleID              
                position = 0       
                Lambda = eye(obj.num_cond)
            end

            x = xy(1) ; y = xy(2) ;
            v_globID = obj.gm.Mesh.Elements(:,eleID) ;
            v = obj.gm.Mesh.Nodes(:,v_globID)' ;
        
            S_loc = det([ones(3,1), v]) ;

            Dxphi = [v(2,2) - v(3,2), ...
                     v(3,2) - v(1,2), ...
                     v(1,2) - v(2,2)] ;
        
            Dyphi = [v(3,1) - v(2,1), ...
                     v(1,1) - v(3,1), ...
                     v(2,1) - v(1,1)] ;
        
            phi   = [det([ ones(3,1) , [ [x,y] ; v(2,:) ; v(3,:) ] ]), ...
                     det([ ones(3,1) , [ v(1,:) ; [x,y] ; v(3,:) ] ]), ...
                     det([ ones(3,1) , [ v(1,:) ; v(2,:) ; [x,y] ] ])] ;

            if position == 0
                coeff = kron([Dxphi ; Dyphi ; phi] ./ S_loc , eye(obj.dep_num));
                res =  [obj.A1(x,y), obj.A2(x,y), obj.A0(x,y)] * coeff ;
            elseif position == 1
                coeff = kron([Dxphi ; Dyphi ; phi] ./ S_loc , eye(obj.dep_num));
                res =  Lambda * [obj.B1(x,y), obj.B2(x,y), obj.B0(x,y)] * coeff;
            end

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
            N_Omega = size(interior_points, 1) ;  % number of interior points
            mesh = obj.gm.Mesh ; 
        
            Cols = cell(1, N_Omega)    ; Rows = cell(1, N_Omega) ;
            B_Omega = cell(1, N_Omega) ; d_Omega = cell(1, N_Omega) ;
        
            % Loop over all interior points and calculate local B_Omega, d_Omega in
            % parallel
        
            hbar = parfor_progressbar_v1(N_Omega,'Local B_{\Omega} calculation, please wait...') ; % progress bar
            parfor i = 1 : N_Omega
                xy = interior_points(i,:) ;
                eleID = obj.findElement(xy, mesh, obj.Omega) ;
        
                if eleID == 0
                   continue
                end
        
                iota = mesh.Elements(:, eleID) ; 
                B_Omega{i} = obj.B_loc(xy, eleID, 0) ;
                d_Omega{i} = obj.f(xy(1),xy(2)) ;
        
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
        
            Cols = reshape(Cols, [], 1) ;
            Rows = reshape(Rows, [], 1) ;
            B_Omega = reshape(B_Omega, [], 1) ;
            B_Omega = sparse(Rows, Cols, B_Omega, maxRow, m*n) ;
            N_Omega_modified = round(maxRow/obj.num_eq) ;
            d_Omega(any(isnan(B_Omega), 2)) = []; B_Omega(any(isnan(B_Omega), 2), :) = [];

            if save_result
                obj.Omega_assemble.B = B_Omega;
                obj.Omega_assemble.d = d_Omega;
                obj.Omega_assemble.N = N_Omega_modified;
            end

        end


        % This function evaluates and assembles global B_Gamma.
        function [B_Gamma, d_Gamma, N_Gamma_modified] = Assemble_B_Gamma(obj, boundary_points, Lambda, save_result)
        
            % boundary_points = [X,Y]
            % Lambda m ny m square(weights) for BC

           arguments 
               obj
               boundary_points
               Lambda
               save_result = true
           end
        
            n = obj.basis_num ;           % number of global basis
            m = obj.dep_num ;             % number of dependent variables
            N_Gamma = size(boundary_points,1) ;   % number of boundary points
            mesh = obj.gm.Mesh ; 
        
            Cols = cell(1, N_Gamma)    ; Rows = cell(1, N_Gamma) ;
            B_Gamma = cell(1, N_Gamma) ; d_Gamma = cell(1, N_Gamma) ;
            
            % Loop over all boundary points and calculate local B_Gamma, d_Gamma in
            % parallel
        
            hbar = parfor_progressbar_v1(N_Gamma,'Local B_{\Gamma} calculation, please wait...') ; % progress bar
            parfor i = 1 : N_Gamma
                xy = boundary_points(i, :);
                eleID = obj.findElement(xy, mesh, obj.Omega);
        
        
                if eleID == 0
                   continue
                end
        
                iota = mesh.Elements(:, eleID) ; 
                B_Gamma{i} = obj.B_loc(xy, eleID, 1, Lambda) ;
                d_Gamma{i} = Lambda*obj.g(xy(1),xy(2)) ;
                
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
            N_Gamma_modified = round(maxRow/obj.num_cond) ;
            d_Gamma(any(isnan(B_Gamma), 2)) = []; B_Gamma(any(isnan(B_Gamma), 2), :) = [];

            if save_result
                obj.Gamma_assemble.B = B_Gamma;
                obj.Gamma_assemble.d = d_Gamma;
                obj.Gamma_assemble.N = N_Gamma_modified;
            end

        end


        % This function evaluate the point-wise numerical value according to the solution coefficient matrix C.
        function u_fit = evaluate_LSFEM_elewise(obj, xy, C, derivative_order)
            % u_fit is a column vector.

            arguments 
                obj
                xy
                C = obj.C
                derivative_order = 0 % or 11 for Dx, 12 for Dy
            end

            if isempty(C)
                C = obj.C;
            end

            try
                mesh = obj.gm.Mesh ;
                eleID = obj.findElement(xy, mesh, obj.Omega) ;
                x = xy(1) ; y = xy(2) ;
                v_globID = mesh.Elements(:,eleID) ;
                v = mesh.Nodes(:,v_globID)' ;

                C_loc = C(:,v_globID) ; S_loc = det([ones(3,1), v]) ;
                v1 = [1, v(1,:)]; v2 = [1, v(2,:)]; v3 = [1, v(3,:)]; vxy  = [1, x, y];
                
                switch derivative_order
                    case 0
                        u_fit = (C_loc(:,1) * det([vxy;v2;v3]) ...
                               + C_loc(:,2) * det([v1;vxy;v3]) ...
                               + C_loc(:,3) * det([v1;v2;vxy])) ./ S_loc;
 
                    case 11
                        vy = v(:,2);
                        u_fit = (C_loc(:,1) * (vy(2)-vy(3)) ...
                               + C_loc(:,2) * (vy(3)-vy(1)) ...
                               + C_loc(:,3) * (vy(1)-vy(2)) ) ./ S_loc;

                    case 12
                        vx = v(:,1);
                        u_fit = (C_loc(:,1) * (vx(3)-vx(2)) ...
                               + C_loc(:,2) * (vx(1)-vx(3)) ...
                               + C_loc(:,3) * (vx(2)-vx(1)) ) ./ S_loc;
                end

            catch 
                u_fit = NaN * C(:,1) ;
            end
        end


        % This function evaluate the vectorized numerical value according to the solution coefficient matrix C.
        function u_fit = evaluate_LSFEM(obj, XY, C, derivative_order, showbar)
            % XY = [ X , Y ]  
            arguments
                obj
                XY
                C = obj.C
                derivative_order = 0
                showbar = true
            end
            
            if isempty(C)
                C = obj.C;
            end

            u_fit = cell(1, size(XY,1)) ;

            if showbar
                hbar = parfor_progressbar_v1(size(XY,1),'Evaluation, please wait...') ;
                parfor i = 1 : size(XY,1)
                    u_fit{i} = obj.evaluate_LSFEM_elewise(XY(i,:), C, derivative_order) ;
                    hbar.iterate(1) ; 
                end
                close(hbar) ;
            else
                parfor i = 1 : size(XY,1)
                    u_fit{i} = obj.evaluate_LSFEM_elewise(XY(i,:), C, derivative_order) ;
                end
            end

           u_fit = cat(2, u_fit{:}) ;
        end


        % This function evaluate the absolute error || u_true(x,y) - u_hat(x,y) ||2 at a series of fixed points.
        function err = abs_error(obj, x, y, C, derivative_order)
            arguments
                obj
                x
                y
                C = obj.C
                derivative_order = 0
            end
            
            if isempty(C)
                C = obj.C;
            end

            [sz1,sz2] = size(x);
            x = reshape(x, [], 1);
            y = reshape(y, [], 1);

            if derivative_order == 0
                err = obj.evaluate_LSFEM([x,y], C, derivative_order, false) - obj.u(x', y');
            elseif derivative_order == 11
                err = obj.evaluate_LSFEM([x,y], C, derivative_order, false) - obj.Dx_u(x', y');
            elseif derivative_order == 12
                err = obj.evaluate_LSFEM([x,y], C, derivative_order, false) - obj.Dy_u(x', y');
            end

            err = sqrt(sum(err .* err, 1))';
            err(isnan(err)) = 0;
            err = reshape(err, sz1, sz2);
        end

        
        % This function approximate the L2 error || u_true - u_hat ||L2.
        function Hk_err = Hk_error(obj, C, type)
            arguments
                obj
                C = obj.C
                type = 0 % Hk error, type = 0 or 1 for k = 0 or 1
            end

            if isempty(C)
                C = obj.C;
            end

            elements = obj.gm.Mesh.Elements;
            nodes = obj.gm.Mesh.Nodes;
            N = size(elements,2);
            val = cell(1,N);

            if type == 0
                abserr2 = @(x,y) obj.abs_error(x,y,C).^2;
            elseif type == 1
                abserr2 = @(x,y) obj.abs_error(x,y,C).^2 + obj.abs_error(x,y,C,11).^2 + obj.abs_error(x,y,C,12).^2;
            end

            parfor eleID = 1:N
                v_globID = elements(:,eleID) ;
                v = nodes(:,v_globID)' ;
                [XY, W] = simplexquad(2,v);
                err = abserr2(XY(:,1), XY(:,2));
                val{eleID} = sum(err .* W);
            end
            Hk_err = sqrt(sum(cell2mat(val)));
        end


        % % This function approximate the L2 error || u_true - u_hat ||L2.
        % function L2_err = L2_error(obj, C)
        %     arguments
        %         obj
        %         C = obj.C
        %     end
        %     elements = obj.gm.Mesh.Elements;
        %     nodes = obj.gm.Mesh.Nodes;
        %     N = size(elements,2);
        %     val = cell(1,N);
        %     abserr = @(x,y) obj.abs_error(x,y,C);
        % 
        %     parfor eleID = 1:N
        %         v_globID = elements(:,eleID) ;
        %         v = nodes(:,v_globID)' ;
        %         [XY, W] = simplexquad(2,v);
        %         err = abserr(XY(:,1), XY(:,2));
        %         val{eleID} = sum((err.^2) .* W);
        %     end
        %     L2_err = sqrt(sum(cell2mat(val)));
        % end


        % This function construct uniform X, Y grids for solution evaluation and visualization.
        function [X, Y, in, XY] = OmegaMeshgrid(obj, grid_size)
            % grid_size = [number of points in x direction, number of points in y direction].
            % in(i,j) returns whether (X(i),Y(j)) is within the domain Omega.
            % XY vectorized version of grids for solution evaluation.
    
            xs = linspace(min(obj.Omega.Vertices(:,1)), max(obj.Omega.Vertices(:,1)), grid_size(1));
            ys = linspace(min(obj.Omega.Vertices(:,2)), max(obj.Omega.Vertices(:,2)), grid_size(2));
            [X,Y] = meshgrid(xs,ys);
            [in, on] = inpolygon(X, Y, obj.Omega.Vertices(:,1), obj.Omega.Vertices(:,2));
            in = in + on;
            in(in==0) = NaN;
            XY = [reshape(X,[],1) , reshape(Y,[],1)];

        end


         % This function reshapes the vectorized c to matrix version C.
         function C = reshape2mat(obj, c)
            C = reshape(c, obj.dep_num, obj.basis_num) ;
         end


         % This function determin whether the point xy is in or on the discretized domain.
         function inon = indomain(obj, xy)
            xy = reshape(xy,2,1) ;
            [in, on] = inpolygon(xy(1), xy(2), obj.Omega.Vertices(:,1), obj.Omega.Vertices(:,2)) ;
            inon = (in | on)*1;
         end


    end



    methods (Static)
        
        % This function samples one point from the boundary of the polygon.
        function sample_coord = sample_polyboundary(boundary_nodes, sample_r)
           % boundary_nodes = [X,Y] : specifies the vertices in a counterclockwise order, 
           %                          with boundary_nodes(1,:) = boundary_nodes(end,:).
           % sample_r : a random float belongs to [0,1].
        
           edge_length = sum((boundary_nodes(1:end-1, :) - boundary_nodes(2:end, :)).^2, 2) ;
           length_prop = cumsum(edge_length) / sum(edge_length) ;
           length_prop = [0;length_prop] ;
                
           sample_edg_ID = min([find(length_prop <= sample_r, 1, 'last'), ...
                               length(edge_length)]) ;
           sample_edg_Nodes = boundary_nodes(sample_edg_ID:sample_edg_ID+1,:) ; % [X Y]
           coord_prop = (sample_r - length_prop(sample_edg_ID))...
                         /(length_prop(sample_edg_ID+1)-length_prop(sample_edg_ID)) ;
           sample_coord = (sample_edg_Nodes'*[1-coord_prop;coord_prop])' ;
         end


        % Interpolate points into the boundary of the polygon.
        function sampled_boundary_nodes = interp_polyboundary(boundary_nodes, para, fixed_Nodes)
           % boundary_nodes = [X,Y] : specifies the vertices in a counterclockwise order, 
           %                          with boundary_nodes(1,:) = boundary_nodes(end,:).

           % para is either a vector rs s.t. rs(i) in [0,1) specifies the relative position of the points ;
           %         or an integer N_gamma specifies the total number of interpolated nodes.

           % fixed_Nodes = [Xs,Ys] specifies the nodes must be included.
           % sampled_boundary_nodes = [Xs, Ys].
        
           arguments 
               boundary_nodes  
               para   
               fixed_Nodes = [] 
           end
        
           if isscalar(para)
               if para-size(fixed_Nodes,1)>= 1
                   rs = linspace(0, 1, para-size(fixed_Nodes,1)) ;
               else
                   return
               end
           end
           
           N = length(rs) ;
           sampled_boundary_nodes = cell(1, length(rs)) ; 
           hbar = parfor_progressbar_v1(N, 'Boundary points sampling, please wait...') ; % progress bar
           parfor i = 1 : N
               sampled_boundary_nodes{i} = Collocation2D.sample_polyboundary(boundary_nodes, rs(i)) ;
               hbar.iterate(1) ;
           end
           close(hbar); % close the progress bar

           sampled_boundary_nodes = [cat(1, sampled_boundary_nodes{:}) ; fixed_Nodes] ;
       end


        % Interpolate points into the boundary of the polygon with holes.
        function sampled_boundary_nodes = interp_multiedge_polyboundary(edge, para, fixed_Nodes)
            % para is a cell of length(edge)

           arguments 
               edge
               para   
               fixed_Nodes = [] 
            end

            num_edg = length(edge) ;

            if ~iscell(para)
                para = repmat({para}, 1, num_edg) ;
            end

            sampled_boundary_nodes = cell(1, num_edg) ;

            for i = 1 : num_edg
                 sampled_boundary_nodes{i} = Collocation2D.interp_polyboundary(edge{i}, para{i}) ;
            end

            sampled_boundary_nodes = [ cat(1, sampled_boundary_nodes{:}) ; fixed_Nodes ] ;

        end


        % This function samples N points from a triangle.
        function sample_coords = sample_from_triangle(XY, N)
        % XY = [X,Y] is a 3*2 matrix specifies the vertices of the triangle.
        % sample_coords = [Xs, Ys] is a N*2 matrix.
    
            sq_r1 = sqrt(rand(N,1)) ; r2 = rand(N,1) ;
            Coes = [1 - sq_r1 , sq_r1.*(1-r2) , sq_r1.*r2] ;
            sample_coords = Coes * XY ;
        end


        % The function additionally samples N points from each element of the triangulated mesh based on the given fixed_Nodes.
        function sampled_interior_nodes = sample_from_triangleMesh(mesh, N, fixed_Nodes)
            % sampled_interior_nodes = [Xs, Ys] is a N*M+n0 by 2 matrix, where M 
            % is the number of the elements, n0 is the muber of nodes specified by fixed_Nodes.
        
            arguments
                mesh   % mesh = gm.Mesh is a mesh object after triangulation.
                N          
                fixed_Nodes = [] % Nodes must be included.
            end
        
            elements = mesh.Elements ; N_ele = size(elements,2) ;
            nodes = mesh.Nodes ;

            sampled_interior_nodes = cell(1, N_ele) ;
            hbar = parfor_progressbar_v1(N_ele, 'Interior points sampling, please wait...') ; % progress bar
            parfor i = 1 : N_ele
                tri = elements(:,i) ;
                triVert = nodes(:,tri)' ;
                sampled_interior_nodes{i} = Collocation2D.sample_from_triangle(triVert, N) ;
                hbar.iterate(1) ;
            end
            close(hbar); % close the progress bar
            sampled_interior_nodes = [cat(1, sampled_interior_nodes{:}) ; fixed_Nodes] ;
         end


         % This function finds the element (ID) where the point xy is located.
        function eleID = findElement(xy, mesh, Omega)

            xy = reshape(xy,2,1) ;
            [in, on] = inpolygon(xy(1), xy(2), Omega.Vertices(:,1), Omega.Vertices(:,2)) ;

            if ~in && ~on
                eleID = 0 ;
                return
            end

            near_node_ID = findNodes(mesh, "nearest", xy) ;
            elemIDs = findElements(mesh, "attached", near_node_ID) ; 
            eleID = elemIDs(1) ;
        
            for id = elemIDs
                tri = mesh.Elements(:,id) ;
                triVert = mesh.Nodes(:,tri)' ;
                [in,on] = inpolygon(xy(1),xy(2),triVert(:,1),triVert(:,2)) ;
                if in || on
                    eleID = id ;
                    return
                end
            end
        end


        % This function reshapes [u1;...;um] to 3D [[u1] ... [um]], [ui] is of size sz1 by sz2. 
        function reshapedArray = reshapeRowsTo3D(A, sz1, sz2)

            [m, ~] = size(A);
            reshapedArray = zeros(sz1, sz2, m);
        
            for i = 1:m
                reshapedArray(:, :, i) = reshape(A(i, :), sz1, sz2);
            end
        end


     end


end
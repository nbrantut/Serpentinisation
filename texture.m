%% 3D representation of the cubes

% 3D representation of the results of the model predicting the bulk
% reaction rate from the micro-mechanics model
clear all
close all
addpath('visualisation/')

%%%%%%%%%%%%%%%%%% parameters for bulk model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = parameters('hsgrain',55e-6,...
    'hsgrainmin',10e-6,...
    'grainsizeprop',2.9,...
    'supcrtfile','SUPCRT/data2_P',...
    'ac0',0.8);

p.sigmainf = 0;

%%%%%%%%%%%%%%%%%% running the bulk model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sol = reaction(p);

%%%%%%%%%%%%%%%%%% extract parameters from the bulk model %%%%%%%%%%%%%%%%%
d = sol.D*p.l*1e6;      % grain diameter
fragm = sol.frag;       % number of fragmentations
fragmu = unique(fragm); % needed for the first loop

%%%%%%%%%%%%%% grain size at which we want to plot stuff %%%%%%%%%%%%%%%%%%
idfig = [45 111 165 210];

%%%%%%%%%%%%%%%%%% position of the cubes as a function of d %%%%%%%%%%%%%%%
X{1} = 0; Y{1} = 0; Z{1} = 0; % initial position of the cube
for i = 2:length(fragmu)
    ind_temp = find(fragm==fragmu(i));
    ind(i) = ind_temp(1); % indice of the first point after fragmentation
    compt_x = 0;
    for ix = 1:size(X{i-1},1)
        compt_x = compt_x + 2;
        compt_y = 0;
        for iy = 1:size(X{i-1},2)
            compt_y = compt_y + 2;
            compt_z = 0;
            for iz = 1:size(X{i-1},3)
                compt_z = compt_z + 2;
                % for each fragmentation, the coordinates of the center of
                % eight cubes are calculated from the position of the
                % pre-existing cube and with side = half of previous cube
                % size
                [X{i}(compt_x-1:compt_x,compt_y-1:compt_y,compt_z-1:compt_z),...
                    Y{i}(compt_x-1:compt_x,compt_y-1:compt_y,compt_z-1:compt_z),...
                    Z{i}(compt_x-1:compt_x,compt_y-1:compt_y,compt_z-1:compt_z)] = ...
                    ndgrid([X{i-1}(ix,iy,iz)-d(ind(i))/2 X{i-1}(ix,iy,iz)+d(ind(i))/2],...
                    [Y{i-1}(ix,iy,iz)-d(ind(i))/2 Y{i-1}(ix,iy,iz)+d(ind(i))/2],...
                    [Z{i-1}(ix,iy,iz)-d(ind(i))/2 Z{i-1}(ix,iy,iz)+d(ind(i))/2]);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%% displaying the cubes as a function of d %%%%%%%%%%%%%%%%

%% definition of the plane intersecting the cubes (plane I)
no_p = [0.2 1  0.7]; % vector normal to plane I
V0   = [5  -1 -3  ]; % one point belonging to plane I
%% coordinates of plane I for plotting on subplot(131)
p_cx = [-p.hsgrain/2 -p.hsgrain/2  p.hsgrain/2 p.hsgrain/2].*10^6;
p_cy = [ p.hsgrain/2 -p.hsgrain/2 -p.hsgrain/2 p.hsgrain/2].*10^6;
p_cz = [(no_p(1)*(p_cx(1)-V0(1))+no_p(2)*(p_cy(1)-V0(2)))/no_p(3)+V0(3) ...
    (no_p(1)*(p_cx(2)-V0(1))+no_p(2)*(p_cy(2)-V0(2)))/no_p(3)+V0(3) ...
    (no_p(1)*(p_cx(3)-V0(1))+no_p(2)*(p_cy(3)-V0(2)))/no_p(3)+V0(3) ...
    (no_p(1)*(p_cx(4)-V0(1))+no_p(2)*(p_cy(4)-V0(2)))/no_p(3)+V0(3)];
%% definition of the matrices of rotation and translation for displaying
%% intersections between the cubes and plane I in 2D (subplot(132))
% line of the intersection of plane I with xy plane
[p_l,no_l,check] = plane_intersect(no_p,V0,[0 0 1],[0 0 0]);
% angle between plane I and xy plane
ang_p = acosd((dot(no_p,[0 0 1]) / norm(no_p)*norm([0 0 1])));
% calculation of the rotation and translation matrices
[R,t] = AxelRot(ang_p,no_l,p_l);

%% plotting the cubes
% defintion of the faces of the cubes
fv.faces = [1 3 2
    1 4 3
    4 5 3
    3 5 6
    7 1 2
    7 2 8
    5 7 6
    6 7 8
    3 6 2
    2 6 8
    5 1 7
    5 4 1]; % cube faces split into triangles defined by the
% vertices number; the order of the vertices determines
% the normal to the faces which is used for Monte Carlo
% loop for plotting everything as d evolves

figure;
k=0;

for i = idfig
    k=k+1;
    
    subplot(2,2,k);
    hold on;
    box on;
    
    clear coord_1x coord_2x coord_1y coord_2y coord_1z coord2_z

    %^ defintion of min and max coordinates for each cube in 3D
    xmin = X{fragm(i)+1}-d(i)/2*ones(size(X{fragm(i)+1}));
    xmax = X{fragm(i)+1}+d(i)/2*ones(size(X{fragm(i)+1}));
    ymin = Y{fragm(i)+1}-d(i)/2*ones(size(Y{fragm(i)+1}));
    ymax = Y{fragm(i)+1}+d(i)/2*ones(size(Y{fragm(i)+1}));
    zmin = Z{fragm(i)+1}-d(i)/2*ones(size(Z{fragm(i)+1}));
    zmax = Z{fragm(i)+1}+d(i)/2*ones(size(Z{fragm(i)+1}));
    
    %% loop on each cube
    for ix = 1:size(X{fragm(i)+1},1)
        for iy = 1:size(X{fragm(i)+1},1)
            for iz = 1:size(X{fragm(i)+1},1)
                
                % vertices of each cube
                fv.vertices = [xmin(ix,iy,iz) ymin(ix,iy,iz) zmax(ix,iy,iz)
                    xmin(ix,iy,iz) ymin(ix,iy,iz) zmin(ix,iy,iz)
                    xmin(ix,iy,iz) ymax(ix,iy,iz) zmin(ix,iy,iz)
                    xmin(ix,iy,iz) ymax(ix,iy,iz) zmax(ix,iy,iz)
                    xmax(ix,iy,iz) ymax(ix,iy,iz) zmax(ix,iy,iz)
                    xmax(ix,iy,iz) ymax(ix,iy,iz) zmin(ix,iy,iz)
                    xmax(ix,iy,iz) ymin(ix,iy,iz) zmax(ix,iy,iz)
                    xmax(ix,iy,iz) ymin(ix,iy,iz) zmin(ix,iy,iz)];

                % coordinates of the cutting edges of the cubes
                fv.edges1 = [fv.vertices(1,:) % first coordinate
                    fv.vertices(1,:)
                    fv.vertices(1,:)
                    fv.vertices(2,:)
                    fv.vertices(2,:)
                    fv.vertices(3,:)
                    fv.vertices(3,:)
                    fv.vertices(4,:)
                    fv.vertices(5,:)
                    fv.vertices(5,:)
                    fv.vertices(6,:)
                    fv.vertices(7,:)];
                fv.edges2 = [fv.vertices(2,:) % second coordinate
                    fv.vertices(4,:)
                    fv.vertices(7,:)
                    fv.vertices(3,:)
                    fv.vertices(8,:)
                    fv.vertices(4,:)
                    fv.vertices(6,:)
                    fv.vertices(5,:)
                    fv.vertices(6,:)
                    fv.vertices(7,:)
                    fv.vertices(8,:)
                    fv.vertices(8,:)];
                
                % coordinates of the intersection of the cutting edges of
                % the cubes with plane I
                [I,check] = plane_line_intersect_vect(no_p,V0,fv.edges1,fv.edges2);
                points_I = I(check==1,:);
                % rotation/translation of the points on the xy plane
                points_I2D = (bsxfun(@plus,R*points_I',t))';
                if length(unique(points_I2D))>1
                    % indices for joining the points
                    ind_ok = convhull(points_I2D(:,1),points_I2D(:,2));
                    % plot the polygons of the intersection
                    fill(points_I2D(ind_ok,1),points_I2D(ind_ok,2),[0.9 0.9 0.9])
                else
                    % plot if not polygons but points only
                    plot(points_I2D(:,1),points_I2D(:,2),'.','Color',[0.9 0.9 0.9])
                end
                % save the initial limit of plane I in 2D for subplot(132)
                if k == 1
                    xlm_save = [min(points_I2D(ind_ok,1)) max(points_I2D(ind_ok,1))];
                    ylm_save = [min(points_I2D(ind_ok,2)) max(points_I2D(ind_ok,2))];
                end
            end
        end
    end
    
    %% second figure with plane I in 2D
    axis equal
    xlim(xlm_save)
    ylim(ylm_save)
    title(['\xi='  num2str(sol.xi(i)*100,3) '%,  \Delta='  num2str(d(i),3) '{\mu}m'])
    % displaying the scale
    len_sc = floor(sum(abs(xlm_save))/5);
    plot([xlm_save(1)+sum(abs(xlm_save))/20;xlm_save(1)+sum(abs(xlm_save))/20+len_sc],...
        [ylm_save(1)+sum(abs(ylm_save))/20;ylm_save(1)+sum(abs(ylm_save))/20],...
        '-k', 'LineWidth', 2)
    text(xlm_save(1)+sum(abs(xlm_save))/20 + len_sc/2,...
        ylm_save(1)+2*sum(abs(ylm_save))/20,...
        [num2str(len_sc) '\mu' 'm'], 'HorizontalAlignment','center')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    set(gca,'Color',[0.7 0.7 0.7]); % color for serpentine in BSE images
    
end

set(gcf, 'Color','w',...
    'InvertHardCopy', 'off');
exportfig('texture','xsize',12,'ysize',13)
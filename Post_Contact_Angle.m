%% Post Processing for Contact angle
% close all

%Finding surface particles.
max_f = sort(curv_a);
distance = sort(sqrt(particles(:,1).^2 + particles(:,2).^2));
max_distance = distance(end);
dist_flag = 0.5*max_distance;
surf_flag = abs(0.6* max_f(1))
n=1;
m=1;

%     figure('Name','Fluid Triangulation')
%     axis([-1 1 0 2])
%     TRI = delaunay(particles(int_fluid,1),particles(int_fluid,2));
%     triplot(TRI,particles(int_fluid,1),particles(int_fluid,2))

figure('Name','Bubble Real Shape')
hold on
% axis([-0.2 0.2 -0.05 .4])
axis([-1 1 -0.1 1.9])
scatter(particles(int_boundary,1),particles(int_boundary,2),[],'filled','MarkerFaceColor',[0.4 0.4 0.4])
% plot(particles(int_boundary,1),particles(int_boundary,2),'.k')
for i = 1:length(int_fluid)
    dist  = sqrt(particles(i,1).^2 + particles(i,2).^2);
    if abs(curv_a(i)) < surf_flag 
        scatter(particles(i,1),particles(i,2),'filled','b');
        surface_indeces(n) = i;
        if particles(i,1)>0
            RHS_SI(m) = i;
            m=m+1;
        end
        n=n+1;
    else
        scatter(particles(i,1),particles(i,2),'filled','g');
    end
end


    x = particles(surface_indeces,1);
    y = particles(surface_indeces,2)
%     figure('Name','Triangulated Surface')
%     TRI = delaunay(x,y);
%     triplot(TRI,x,y)
%     
    figure
    hold on
    scatter(particles(int_boundary,1),particles(int_boundary,2),'k','filled','S')
%     y(end) = NaN;
%     c = PF_sum(surface_indeces);
%     patch(x,y,c,'EdgeColor','interp','Marker','o','MarkerFaceColor','flat');
    
    
    
%     figure('Name','Slopes to root')
%     hold on
%     plot(particles(int_boundary,1),particles(int_boundary,2),'.k')
    x = particles(RHS_SI,1);
    y = particles(RHS_SI,2);
    length(y)
    plot(x,y,'o')

%     figure
%     y(end) = NaN;
%     c = PF_sum(RHS_SI);
%     patch(x,y,c,'EdgeColor','interp','Marker','o','MarkerFaceColor','flat');
    
    % Sort by Ascending Y-coord
    surf_particles = particles(RHS_SI,:);
    [sort_y,sort_y_I] = sort(y);
    length(sort_y_I)
    root_x = surf_particles(sort_y_I(1),1)
    root_y = surf_particles(sort_y_I(1),2)
    plot(root_x,root_y,'ro')
    
    % Sort by Ascending x-coord
    surf_particles = particles(RHS_SI,:);
    [sort_x,sort_x_I] = sort(x);
    length(sort_x_I)
    root_x = surf_particles(sort_x_I(end),1)
    root_y = surf_particles(sort_x_I(end),2)
    plot(root_x,root_y,'mo')
%     axis([-1 1 0 2])
    
    for j = 2:length(RHS_SI)

        leaf_x = x(sort_y_I(j));
        leaf_y = y(sort_y_I(j));
        diff = -[root_x ; root_y] + [leaf_x ; leaf_y];
        slope(j) = (leaf_y - root_y) / (leaf_x - root_x);
        line([root_x leaf_x], [root_y  leaf_y])
        
    end
    slope_avg = (sum(slope) - max(slope) + min(slope))/length(slope);
    Contact_angle = 180-abs(atand(slope_avg))
    
%     figure('Name','Traingulated RHS Surface')
%     TRI = delaunay(x,y);
%     triplot(TRI,x,y)
%     

% finding boundary
    boundy=max(particles(int_boundary,2));
    boundy=(boundy+0.5*dx)*ones(1,100);
    boundx=linspace(min(particles(int_boundary,1)),max(particles(int_boundary,1)),100);
    
figure('Name','Color function gradient')
hold on
grid on
box on
grid minor
colormap jet
c=curv_a;
scatter(particles(1:length(curv_a),1),particles(1:length(curv_a),2),[],c,'filled');
scatter(particles(int_boundary,1),particles(int_boundary,2),[],'filled','MarkerFaceColor',[0.4 0.4 0.4])
plot(boundx,boundy,'k','LineWidth',2);
axis([-0.2 0.2 -0.05 0.4])
hcb=colorbar;
title(hcb,'$color gradient$','FontSize',14,'Interpreter','Latex');

% xxx=-0.05*ones(1,100);
% yyy=linspace(0,0.1,100);
% plot(xxx,yyy,'k','LineWidth',2);

% ref_theta=30;
% refx=-linspace(0,0.1,1000);
% refy=-refx*sind(ref_theta);
% plot(root_x+0.5*dx+refx,refy+0.5*dx,'k','LineWidth',2);
        
        
        
        
        
        
        
        
        
        
        
        
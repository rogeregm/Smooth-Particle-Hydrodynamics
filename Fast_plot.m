function [] = Fast_plot(particles1,int_fluid,p_ref)

% figure
%  hold on
% axis([-0.2 0.2 -0.1 0.4])
% for i=1:max(size(particles))
%    if particles(i,3) == 1
%        plot(particles(i,1),particles(i,2),'r.')
%    elseif particles(i,3) == 2
%        plot(particles(i,1),particles(i,2),'g.')
%    elseif particles(i,3) == 3
%        plot(particles(i,1),particles(i,2),'m.')
%    end
% end

avg_vel_FP = sum(sqrt(particles1(int_fluid,6).^2 + particles1(int_fluid,7).^2))/length(int_fluid);
% figure('Name','Velocity plot')
% axis([-1 1 -0.2 2])
% quiver(particles(:,1),particles(:,2),particles(:,6),particles(:,7));
% text(0,-0.1,['Average Velocity: ',num2str(avg_vel_FP)]);

figure('Name','Pressure Distribution')
hold on
box on
set(gcf, 'Position', 540*[0.2, 0.2, 1.25, 0.75]);
axis([-0.56 0.56 -0.07 1.12])

%     load('save_var_crown2.mat')
%     scatter(particles(int_boundary,1),particles(int_boundary,2),'filled','MarkerFaceColor',[0.4 0.4 0.4]);
    
colormap jet
c=particles1(int_fluid,8)/p_ref(2);
% caxis([0 1.2])
wall_size=(length(int_fluid)+1:length(particles1));
scatter(particles1(int_fluid,1),particles1(int_fluid,2),60,c,'filled')
scatter(particles1(wall_size,1),particles1(wall_size,2),60,'filled','MarkerFaceColor',[0.4 0.4 0.4])
hcb=colorbar
title(hcb,'$p/p_{ref}$','FontSize',14,'Interpreter','Latex')


figure('Name','Velocity Distribution')
hold on
box on
set(gcf, 'Position', 540*[0.2, 0.2, 1.25, 0.75]);
axis([-0.56 0.56 -0.07 1.12])

%     load('save_var_crown2.mat')
%     scatter(particles(int_boundary,1),particles(int_boundary,2),'filled','MarkerFaceColor',[0.4 0.4 0.4]);

    
colormap jet
c=sqrt(particles1(int_fluid,6).^2+particles1(int_fluid,7).^2);
% caxis([0 2.2789])
wall_size=(length(int_fluid)+1:length(particles1));
scatter(particles1(int_fluid,1),particles1(int_fluid,2),60,c,'filled')
scatter(particles1(wall_size,1),particles1(wall_size,2),60,'filled','MarkerFaceColor',[0.4 0.4 0.4])
hcb=colorbar;
title(hcb,'$v$','FontSize',14,'Interpreter','Latex')
end
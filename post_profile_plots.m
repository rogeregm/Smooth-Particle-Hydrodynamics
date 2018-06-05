%% profile plots

Height_fs = max(particles(int_fluid,2))+h/2;
y = 0: Height_fs/20 : Height_fs;
p_ana = rho_0(2) * abs(g(2)) * (Height_fs - y)/p_ref(2);

figure
fz = 16;
width = 1;
hold on
box on
grid on
set(gcf, 'Position', (540*width + 20)*[0.2, 0.2, 1.0, 0.40]);
plot(y,p_ana,'color',[0 0.5 1],'linewidth',3)
% for n = 100:104
% scatter(particles(int_fluid,2,n),particles(int_fluid,8,n)/p_ref(2),'.k')
scatter(particles(int_fluid,2),particles(int_fluid,8)/(p_ref(2)),'.k')
% end
leg=legend('Analytical','Simulation')
legend('boxoff')
set(leg,'FontSize',14,'Interpreter','Latex')
xlabel('Height', 'fontsize',14,'Interpreter','Latex')
ylabel('$p/p_{ref}$','fontsize',14,'Interpreter','Latex')
set(gca, 'fontsize',14)


%% test 1 and 2
% pressure profile
%{
kernel = 1;
d = 2;
dx = 0.02;
h = dx;
G = 1 * [ 0, -9.81];
rho_0(2) = 1;
H = max(particles(int_fluid,2)) - min(particles(int_fluid,2)) + dx; % height of fluid
v_ref = sqrt(abs(G(2)) * H );
t_ref = H / v_ref; 
p_ref = rho_0 * abs(G(2)) * H;
%}

average = 0; % 1 yes
if average == 1
    time = 200:500;
    particles(int_fluid,:) = mean(save_pos(int_fluid,:,time),3);
end
H_max = max(particles(int_fluid,2))+h/2;
H_min = min(particles(int_fluid,2))-h/2;
y =  H_min : abs(H_max - H_min)/100 : H_max ;
x = 0;
press(1 : size(y,2)) = 0;

% cutoff radius
if kernel == 1
    r_c = 3 * h;
elseif kernel == 2
    r_c = 2 * h;
end
% range-search
Profile(:,2) = [x,y];
NS = particles(int_fluid,1:2);
[idx] = rangesearch(NS,Profile,r_c);

for n = 1 : size(y,2)
    
    sum_rhoW = 0;
    sum_pW = 0;
    sum_W = 0;
    
    for m = idx{n}
        
        %distance between particles
        drx = x - particles(m,1);
        dry = y(n) - particles(m,2);
        rad = sqrt(drx^2 + dry^2);

        % nondimensional distance
        q = rad / h;
        W = kernel_fct(kernel, d, h, q);
        sum_rhoW = sum_rhoW + particles(m,5) * W;
        sum_pW = sum_pW + particles(m,8) * W +PF_sum(m)*W;
        sum_W = sum_W + W;        
    end
    press1(n) = sum_pW / sum_W;
    rho2(n) = sum_rhoW / sum_W;
end
press2 = p_0(2) * ( (rho2/rho_0(2)).^gamma(2) - 1 ) + Xi(2);

figure
hold on
grid on
grid minor
set(gcf, 'Position', (540 + 20)*[0.2, 0.2, 1.0, 0.40]);

plot(y,(press1)/(p_ref(2) + Xi(2)),'ok')
% plot(y,(press2)/(p_ref(2) + Xi(2)),'.k')
plot(y,rho_0(2) * abs(G(2)) * (H-y)/(p_ref(2) + Xi(2)), 'color',[0 0.5 1],'LineWidth',2)

figure
plot(y,press1-press2)

figure('Name','Pressure Distribution')
hold on
box on
set(gcf, 'Position', 540*[0.2, 0.2, 1.25, 0.75]);
axis([-0.56 0.56 -0.07 1.12])

%     load('save_var_crown2.mat')
%     scatter(particles(int_boundary,1),particles(int_boundary,2),'filled','MarkerFaceColor',[0.4 0.4 0.4]);
    
colormap jet
pdiff=((particles(int_fluid,8)-rho_0(2) * abs(G(2)) * (H-particles(int_fluid,2)))/p_ref(2));
c=pdiff;
caxis([-0.08 0.08])
wall_size=(length(int_fluid)+1:length(particles));
scatter(particles(int_fluid,1),particles(int_fluid,2),60,c,'filled')
scatter(particles(wall_size,1),particles(wall_size,2),60,'filled','MarkerFaceColor',[0.4 0.4 0.4])
hcb=colorbar

title(hcb,'$p/p_{ref}$','FontSize',14,'Interpreter','Latex')


%% scatter profile and time

snaps = 860:1000;
pos = 51:100:4451;
X(:,snaps-min(snaps)+1) = save_pos(pos,1,snaps);
Y(:,snaps-min(snaps)+1) = save_pos(pos,2,snaps);
P(:,snaps-min(snaps)+1) = save_pos(pos,8,snaps);

figure
hold on
for n = 1:45
scatter(Y(n,:),P(n,:))
end




%% Crown Post Processing

%% Tracking crown
close all
load('save_var_crown2.mat')
addpath('/home/reg216/homedir/Project/MatlabSPH_Clean/MatlabSPH_Clean/particle_results_crown2')

frames=31:80;%frames
nn=1;
for n=frames
    if n<10
        filename=['save_pos00',num2str(n)];
    elseif n<100
        filename=['save_pos0',num2str(n)];
    else
        filename=['save_pos',num2str(n)];
    end 

    % load fluid
    clear particles
    load([filename,'.mat'])
    particles=eval(filename);

    [r(n) ,I(n)]=max(particles(:,2));

    if particles(I(n),1)>0
        xcoords(nn)=particles(I(n),1);
        nn=nn+1;
    end
    
    [B,I]=sort(particles(:,2),'descend');
    
    if particles(I(1),1)>0
        xcoord(n)=particles(I(1),1);
    elseif particles(I(1),1)<0
        if particles(I(2),1)>0
            xcoord(n)=particles(I(2),1);
        elseif particles(I(2),1)<0 && particles(I(3),1)>0
            xcoord(n)=particles(I(3),1);
        else
            fprintf('Not in top 3');
        end
    end
            
    cons_mass(n)=sum(particles(:,4));
    clear filename
end

figure('Name','x')
grid on
grid minor
box on
plot(2.2637*10^-4*(1:length(xcoords))./sqrt((1000*0.2^3)/0.019067),xcoords./0.2,'.')

figure('Name','y')
hold on
grid on
grid minor
plot(log(save_pos_t(frames)*5/0.2),log(r(frames)./0.2),'o','MarkerFaceColor',[0 0 1])
p1 = polyfit(log(save_pos_t(frames)*5/0.2),log(r(frames)./0.2),1)
y1=polyval(p1,log(save_pos_t(frames)*5/0.2));
plot(log(save_pos_t(frames)*5/0.2),y1,'--g','linewidth',2)
square_root_time_1=sqrt(save_pos_t(frames(end)))

%% Tracking base

addpath('/home/reg216/homedir/Project/MatlabSPH_Clean/MatlabSPH_Clean/particle_results_crown2')

frames=31:80;%frames


%test area : base definition
profile(1:50,1)=linspace(0,0.5,50);
profile(:,2)=0.4;

for n=frames
    %load file
    if n<10
        filename=['save_pos00',num2str(n)];
    elseif n<100
        filename=['save_pos0',num2str(n)];
    else
        filename=['save_pos',num2str(n)];
    end 

    % load fluid
    clear particles
    load([filename,'.mat'])
    particles=eval(filename);
    
    % find particles within one particle volume of material point (base)
    idx=rangesearch(particles(:,1:2),profile,h);
    fprintf(['frame number ',num2str(n),'\n'])
    
    nn=1;
    for nn=1:length(profile(:,1))
        if isempty(idx{nn})==0
            rbase=particles(idx{nn}(1),1);
            fprintf(['broke at',num2str(nn),'\n']);
            break
        end
    end
    nn
    base_radius(n,1)=n;
    base_radius(n,2)=rbase;
    clear filename
end

% base_radius=unique(base_radius,'stable');
% plot(1:length(rbase),log(rbase/0.2),'.')
% plot(1:length(base_radius),log(base_radius),'.')
% time=2.2637*10^-4.*(31:120)*5/0.2;
% figure('Name','Base Radius')
% hold on
% grid on
% grid minor
% plot((time),(base_radius(31:120,2)/0.2),'.')

%% ploting
figure('Name','Base Radius')
hold on
grid on 
grid minor
box on

%Base
plot(log(save_pos_t(frames)*5/0.2),log(base_radius(frames,2)./0.2),'o','MarkerFaceColor',[0 0 1])
%tip
plot(log(save_pos_t(frames)*5/0.2),log(r(frames)./0.2),'o','MarkerFaceColor',[1 0 1])
%polyfit base
p = polyfit(log(save_pos_t(frames)*5/0.2),log(base_radius(frames,2)./0.2)',1)
y11=polyval(p,log(save_pos_t(frames)*5/0.2));
plot(log(save_pos_t(frames)*5/0.2),y11,'--g','linewidth',2)

%Polyfit Tip
p1 = polyfit(log(save_pos_t(frames)*5/0.2),log(r(frames)./0.2),1)
y1=polyval(p1,log(save_pos_t(frames)*5/0.2));
plot(log(save_pos_t(frames)*5/0.2),y1,'--k','linewidth',2)
square_root_time_1=sqrt(save_pos_t(frames(end)))

xlabel('$\log( \hat {t})$','Interpreter','Latex','Fontsize',14)
ylabel('$\log(r/D)$','Interpreter','Latex','Fontsize',14)
lg=legend('$r_b$','$r_c$')
legend('boxoff')
set(lg,'Interpreter','Latex','Fontsize',14,'Location','northoutside','Orientation','horizontal')

square_root_time=sqrt((save_pos_t(frames(end))))%-save_pos_t(frames(1))))



average_slope=(p(1)+p1(1))*0.5




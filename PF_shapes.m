%% PairWise force graph

range = [0:0.001:1];
h=1;
k=1;

cosine_F = cos(((1.5*pi) / (k*h)).*range);

eps=1/3.5;
eps_0=eps/2;
psi_eps=exp(-range.^2/(2*eps^2));
psi_eps_0 = exp(-range.^2/(2*eps_0^2));

A_2 = (eps/eps_0)^2;
PF_2 = (-A_2.*psi_eps_0 + psi_eps);

A_3 = (eps/eps_0)^3;
PF_3 = -range.*(-A_3.*psi_eps_0 + psi_eps);


for i= 1:length(range)
    
    B_T(i) = kernel_fct(1, 2, 3*0.05, range(i)/(3*0.05));
    
    if 2*range(i) > k*h && range(i) < k*h
        Cohesion_F(i)=(h-range(i))^3*range(i)^3;
    elseif 2*range(i) <= k*h && range(i) >= 0
        Cohesion_F(i)=2*(h-range(i))^3*range(i)^3 - (h^6)/64;
    else
        Cohesion_F(i)= 0;
    end

    if 2*range(i) > k*h && range(i)< k*h
        Adhesion_F(i) = (-(4*range(i)^2)/h + 6*range(i) - 2*h)^(1/4);
    else
        Adhesion_F(i) = 0; 
    end
end
B_T_norm = -B_T/max(B_T);
Cohesion_F = Cohesion_F.*32/(pi*(k*h)^9);
Cohesion_F_norm = Cohesion_F./max(Cohesion_F);
Adhesion_F = 0.007/(h^3.25).*Adhesion_F;
Adhesion_F_norm = Adhesion_F./max(Adhesion_F);

figure
hold on
box on
plot(range,-cosine_F,'k')
plot(range,1/3*PF_2,'g')
plot(range,-PF_3,'r')
plot(range,Cohesion_F_norm,'b')
plot(range,B_T_norm,'c')
plot(range,zeros(1,length(range)),'--');
title('Inter-Particle Force models')
legend( ['Cosine Force: ' num2str(max(cosine_F))],['PF 2: ' num2str(max(PF_2))],['PF 3: ' num2str(max(PF_3))], ['Cohesion Force: ' num2str(max(Cohesion_F))], ['Becker and Teschner: ' num2str(max(B_T))] ,'Location','southeast');
legend('boxoff')
xlabel('Range')
ylabel('Normalized force')
text(0.025,0.1,'\bf Attractive','Fontsize',12)
text(0.025,-0.1,'\bf Repulsive','Fontsize',12)
hold off

figure
plot(range,Cohesion_F_norm,'b')
hold on
plot(range,Adhesion_F_norm,'r')
plot(range,zeros(1,length(range)),'--');
% title('Cohesion and Adhesion Forces')
ylabel('Normalized Force','Fontsize',14,'Interpreter','latex')
xlabel('Range','Fontsize',14,'Interpreter','latex')
leg=legend(['Cohesion Force: ' num2str(max(Cohesion_F))],['Adhesion Force: ' num2str(max(Adhesion_F))],'Location','southeast')
legend('boxoff')
set(leg,'FontSize',14,'Interpreter','latex')
text(0.025,0.1,'\bf Attractive','Fontsize',14,'Interpreter','latex')
text(0.025,-0.1,'\bf Repulsive','Fontsize',14,'Interpreter','latex')



%%
%This code simulates the FDM signal using the approach in:
%Tendler, Benjamin C., and Richard Bowtell. "Frequency difference mapping applied to the corpus callosum at 7T." Magnetic resonance in medicine 81.5 (2019): 3017-3031.
%This code recreates Figure 1 in the paper
%%
%For the simulation, compartment parameters taken from:
%Sati, Pascal, et al. "Micro-compartment specific T2? relaxation in the brain." Neuroimage 77 (2013): 268-278. (Table 1 - SCC)
%%
%Set parameters
%TEs
TE=[2.4:2.4:48]/1000;
%Compartment 1
freq_1=31.8;
T2_1=1/164.2;
A_1=0.126;
%Compartment 2
freq_2=2;
T2_2=1/39.5;
A_2=0.501;
%Compartment 3
freq_3=-4.1;
T2_3=1/24.1;
A_3=0.373;
%Time independent phase offset and average frequency offset
phi_0=0.5;
Omega=50;
%%
%Simulate 3-compartment signal
S=exp(1i*phi_0).*exp(1i*2*pi*Omega.*TE).*(A_1.*exp(-TE./T2_1).*exp(1i*2*pi*freq_1.*TE)+A_2.*exp(-TE./T2_2).*exp(1i*2*pi*freq_2.*TE)+A_3.*exp(-TE./T2_3).*exp(1i*2*pi*freq_3.*TE));
%Convert to phase
pha_s=angle(S);
%%
%Perform FDM steps

%Eq. [3]
pha_s_1=angle(exp(1i*pha_s(:,2:end))./repmat(exp(1i*pha_s(:,1)),[1 length(TE)-1]));

%Eq. [4]
for k=1:size(pha_s_1,2)
    pha_s_2(:,k)=angle(exp(1i*pha_s_1(:,k))./(exp(1i*k*pha_s_1(:,1))));
end
   
%Eq. [5]
for k=2:size(pha_s_2,2)
    FDM(:,k)=pha_s_2(:,k)/(2*pi*(k-1)*(TE(2)-TE(1)));    
end       
%%
%Perform Noise calculation
%Define mean R2* as equal to 1/30ms
mean_R2=1/0.03;
%Calculate echo spacing
delta_TE=TE(2)-TE(1);
%Define SNR
SNR_1=300;

%Eq. [7]
for k=3:size(pha_s,2)
    noise(k)=1/(2.*pi.*delta_TE.*SNR_1)*((exp((k-1)*mean_R2*delta_TE)/(k-2))^2+1+((k-1)*exp(mean_R2*delta_TE)/(k-2))^2)^0.5;
end
%%
%Plot figures
%Phase
figure('units','centimeters','position',[16.8*2 12.6*2 16.8*2 12.6*2]);
subplot(2,2,1)
plot(TE*1000,pha_s,'-o','LineWidth',2,'MarkerFaceColor',[0 0.4470 0.7410],'MarkerSize',4); 
axis([-0 TE(end)*1000+2 -pi pi]);
xlabel('TE (ms)','interpreter','latex') % x-axis label
ylabel('$\phi$ (rad)','interpreter','latex') % y-axis label
set(gca,'YTick',([-pi 0 pi]))
set(gca,'yticklabels',{'$-\pi$','$0$','$\pi$'},'TickLabelinterpreter','latex')
set(gca,'TickLabelinterpreter','latex')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
%Phase after first division
subplot(2,2,2);
plot(TE(2:end)*1000,pha_s_1,'-o','LineWidth',2,'MarkerFaceColor',[0 0.4470 0.7410],'MarkerSize',4); 
axis([-0 TE(end)*1000+2 -pi pi]);
xlabel('TE (ms)','interpreter','latex') % x-axis label
ylabel('$\phi''$ (rad)','interpreter','latex') % y-axis label
set(gca,'YTick',([-pi 0 pi]))
set(gca,'yticklabels',{'$-\pi$','$0$','$\pi$'},'TickLabelinterpreter','latex')
set(gca,'TickLabelinterpreter','latex')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
%Phase after second division
subplot(2,2,3);
plot(TE(3:end)*1000,pha_s_2(2:end),'-o','LineWidth',2,'MarkerFaceColor',[0 0.4470 0.7410],'MarkerSize',4); 
axis([-0 TE(end)*1000+2 -0.9 0]);
xlabel('TE (ms)','interpreter','latex') % x-axis label
ylabel('$\phi''''$ (rad)','interpreter','latex') % y-axis label
set(gca,'TickLabelinterpreter','latex')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
%FDM signal
subplot(2,2,4)
plot(TE(3:end)*1000,FDM(2:end),'-o','LineWidth',2,'MarkerFaceColor',[0 0.4470 0.7410],'MarkerSize',4); 
axis([-0 TE(end)*1000+2 -3 -1.5]);
xlabel('TE (ms)','interpreter','latex') % x-axis label
ylabel('Change in frequency (Hz)','interpreter','latex') % y-axis label
set(gca,'TickLabelinterpreter','latex')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
%Noise simulation
figure('units','centimeters','position',[16.8 12.6 16.8 12.6]);
plot(TE(3:end)*1000,noise(3:end),'-o','LineWidth',2,'MarkerFaceColor',[0 0.4470 0.7410],'MarkerSize',4); 
axis([-0 TE(end)*1000+2 0.3 0.6]);
xlabel('TE (ms)','interpreter','latex') % x-axis label
ylabel('Noise (Hz)','interpreter','latex') % y-axis label
set(gca,'YTick',([0.3 0.4 0.5 0.6]))
set(gca,'TickLabelinterpreter','latex')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
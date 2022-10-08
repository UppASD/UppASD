%Plots the adiabatic magnon spectra W(Q)
%Reads the ams-file

% Clear variables
clear
load ams.SDMTRLNO.out
spindisp_SDMTRLNO = ams_SDMTRLNO;
load qpoints.out

[n1,n2]=size(spindisp_SDMTRLNO(:,:))
figure(12); clf
plot(qpoints(:,1),spindisp_SDMTRLNO(:,2:n2-1),'.-','LineWidth',2)
xlabel('Wavevector Q','FontSize',16)
ylabel('w(Q) meV','FontSize',16)
box on
set(gca,'FontSize',16)
set(gca,'LineWidth',2)
axis tight
print -depsc2 spindisp

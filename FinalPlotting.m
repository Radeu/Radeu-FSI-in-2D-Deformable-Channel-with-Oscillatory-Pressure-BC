clear all;
close all;
format long;
FoldersDir = '/home/tools/a/urade/Downloads/';
FoldersList = {'FSIGm0pt5Wo1'};
Lengths = [1.0256e-3,1.4504e-3,1.7763e-3];
%control_new ="GmSameWoVary";
%control_new = "WoSameGmVary";
%control_new = "P0";
%control_new = "UxUy";
%control_new = "AvgUx";
%control_new = "Pressure";
%control_new = "Velocity";
%control_new = "Mesh";
%control_new = "AvgUy";
%control_new = "U0";
control_new = "U0aP0aAppendix";

if control_new=="GmSameWoVary"
    for i=1:length(FoldersList)
        CurrentFolder = FoldersList{i};
        ModifiedFolder = regexprep(CurrentFolder, 'pt', '.');
        numbers = regexp(ModifiedFolder, '\d+(\.\d+)?', 'match');
        numeric_values = str2double(numbers);
        Gamma = numeric_values(1);
        Wo = numeric_values(2);
        if Gamma==0.5
            L = Lengths(1);
        elseif Gamma==1.0
            L = Lengths(2);
        elseif Gamma==1.5
            L = Lengths(3);
        end
    
        p0 = 7584.3;
         h0 = 4.5e-5;
         b = 4.5e-5; 
         E = 5.9e+5;
         nu_s = 0.45;
         omega = 17951.94;%5714*pi;
         epsilon = h0/L;
         G=E/(2*(1+nu_s));
         K=E/(3*(1-2*nu_s));
         lamda = K-2/3*G;
         kappa_param = b/((lamda+2*G));
         beta = kappa_param*p0/h0;
         defscale = beta*h0;
        C_W = (1-2*nu_s)/(2*(1-nu_s))*b/G;
        C_I = (2*nu_s)*(nu_s-1/4)/(3*(1-nu_s)^2)*b^3/G;
        theta = C_I/(C_W*L^2);
    
        %% 
        pmean_files = dir(append(FoldersDir,CurrentFolder,'/pmean*.csv'));
        Z = linspace(0,1,500);
        AllData = [];
    
        for iii=1:350
            T = readtable(append(FoldersDir,CurrentFolder,'/',pmean_files(iii).name));
            AllData = [AllData,T{:,:}];  
        end
        TempAvgP = mean(AllData,2);
        fig1=figure(1);
        colorarray_sim=viridis(length(FoldersList));
        count_sim = i;
        set(fig1,'DefaultAxesFontName','Times New Roman')
        set(fig1,'DefaultAxesFontSize',18)
        h(i)=plot(Z(1:15:end),TempAvgP(1:15:end)/p0/(beta),'Square','Color',colorarray_sim(count_sim,:), 'MarkerSize',7,'MarkerFaceColor',colorarray_sim(count_sim,:));
        set(gca,'defaulttextinterpreter','latex')
        set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
        set(gca,'box','on');
        xlabel('$Z$');
        ylabel('$\langle P_1 \rangle$');
    
        hold on;
    
        %%
        syms kappa
        syms Womer Gm Th Y Z_ZR Tm real positive
        
        control = "pressure";
        Q1at0 = 0; 0.0019;
        
        % leading-order pressure
        if control == "pressure"
            P0a = sinh(kappa*(1-Z_ZR))/sinh(kappa);
        elseif control == "flow"
            P0a = sinh(kappa*(1-Z_ZR))/(kappa*cosh(kappa));
        else
            error(' ERROR: incorrect control type, flow or pressure');
        end
        P0 = P0a*exp(1i*Tm);
        
        % leading-order axial velocity
        Vz0a = 1/(1i*Womer^2)*(1 - cos(1i^(3/2)*(1-2*Y)*Womer/2)/cos(1i^(3/2)*Womer/2))*(-diff(P0a,Z_ZR));
        Vz0 = Vz0a*exp(1i*Tm);
        
        % leading-order vertical velocity, found from COM
        syms f(Y)
        ode = diff(f,Y) + diff(Vz0a,Z_ZR) == 0;
        cond = f(0) == 0;
        Vy0a = dsolve(ode,cond);
        Vy0 = Vy0a*exp(1i*Tm);
        
        % pre-compute the advective streaming terms
        inertia_avg_term1 = simplify(1/2*real(conj(Vy0)*diff(Vz0,Y)));
        inertia_avg_term2 = simplify(1/2*real(conj(Vz0)*diff(Vz0,Z_ZR)));
        
        % now define the auxilliary function -- if we defined them earlier it
        % makes the symbolic calculation harder
        f0 = (-1j*Womer - 2*1j^(3/2)*tan(1j^(3/2)*Womer/2))/Womer^3;
        kappa = sqrt(1j*Gm/(f0 + Th*1j*Gm));
        % substitute kappa and f0 into the inertia_avg_term
        inertia_avg_term1 = subs(inertia_avg_term1);
        inertia_avg_term2 = subs(inertia_avg_term2);  
        
        % numerically  
        opts = bvpset('RelTol',1e-6,'AbsTol',1e-6,'Stats','on');
        ymesh = linspace(0,1,100);
        solinit = bvpinit(ymesh, @(x)([x*(1-x); 1-2*x]));
        bcfun = @(ya,yb)([ya(1); yb(1)]);
        
        % need to set these for definiteness
        Gm = Gamma;
        Womer = Wo;
        Th = theta;
        
        % now solve for int <V_{z,1}^adv>(Y,Z) dY at each Z in zmesh
        zmesh = linspace(0,1,50);
        yint_Vz1_adv = 0*zmesh;
        for ii = 1:length(zmesh)
            % get current Z value and evaluate inertial_avg_termws
            Z_ZR = zmesh(ii);
            inertia_avg_term1_2 = subs(Womer^2/Gm*inertia_avg_term1);
            inertia_avg_term2_2 = subs(Womer^2/Gm*inertia_avg_term2);
            
            % define ODE RHS, taking into account weird behavior when 
            % argument of matlabFunction is zero
            if eval(inertia_avg_term1_2 + inertia_avg_term2_2 == 0)
                oderhs = @(Y)(0*Y);
            else
                oderhs = matlabFunction(inertia_avg_term1_2 + inertia_avg_term2_2);
            end
            bvpfun = @(Y,y)([y(2); oderhs(Y)]); 
            
            % find solution to BVP for <V_{z,1}^adv>(Y) at this Z location
            sol = bvp4c(bvpfun,bcfun,solinit,opts);
            
            % integrate from 0 to 1 dY and save it
            yint_Vz1_adv(ii) = trapz(sol.x,sol.y(1,:));
        
        end
        
        % now evaluate the streaming flow rate <Q_1> = const.
        syms Z_ZR % first clear Z value
        % this is the displacement based on the combined foundation model
        U0 = P0 - Th*diff(P0,Z_ZR,2);
        % this is the "slip" velocity, <U_0 dV_{z,0}/dY> at Y = 1
        U0dVz0dY = matlabFunction(eval(subs(1/2*real(conj(U0)*subs(diff(Vz0,Y),Y,1)))));
        temp = -1/2*U0dVz0dY(zmesh) + yint_Vz1_adv;
        
        % now evaluate the streaming pressure <P_1>(Z) 
        if control == "pressure"
            % such that <P_1>(0) = <P_1>(1) = 0
            Q1 = trapz(zmesh,temp);
            P1 = 12*cumtrapz(zmesh,temp) - 12*zmesh*Q1;
        elseif control == "flow"
            % such that <P_1>(1) = 0 and <Q_1> is given, maybe = 0
            fzmesh = flip(zmesh);
            ftemp = flip(temp);
            %Q1 = temp(end)/10*ones(size(fzmesh)); 
            Q1 = Q1at0*ones(size(fzmesh));
            %load TimeAvgQ1.mat
            %Q1 = interp1(linspace(0,1,100),norm_TempAvgQ1,fzmesh);    
            P1 = 0*ftemp;    
            % integral from 1 to Z
            for ii=1:length(zmesh)-1
                idx = length(zmesh)+1-ii;
                P1(ii) = 12*trapz(fzmesh(1:idx),ftemp(1:idx)) ...
                         + 12*(1-zmesh(ii))*Q1(idx);
            end
            P1(end) = 0;
        
        else
            error(' ERROR: incorrect control type, flow or pressure');
        end
        
        
        % % make a nice plot
        colorarray_theory=viridis(length(FoldersList));
        count_theory = i;
        set(fig1,'DefaultAxesFontName','Times New Roman')
        set(fig1,'DefaultAxesFontSize',18)
        plot(zmesh,P1,'Color',colorarray_theory(count_theory,:),'LineWidth',1.5);
        set(gca,'defaulttextinterpreter','latex')
        set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
        set(gca,'box','on');
        xlabel('$Z$');
        ylabel('$\langle P_1 \rangle$');
        hold on;
        Wos(i) = Wo;
    end
    ylim([-0.1,0.7]);
    legendEntries = arrayfun(@(wo) ['$\mathrm{Wo} = ' num2str(wo) '$'], Wos, 'UniformOutput', false);
    leg=legend(h,legendEntries, 'Interpreter', 'latex');
    set(leg, 'Location', 'south');
    leg.Position(1) = 0.5 - leg.Position(3)/2; 
    leg.Position(2) = 0.22; 
    GammaStr = strrep(sprintf('%.1f', Gamma), '.', 'pt');
    exportgraphics(figure(1), sprintf('/home/tools/a/urade/Downloads/Feb24/Gm_%s_Wos.pdf', GammaStr), ...
         'ContentType', 'vector');    

elseif control_new=="WoSameGmVary"
    for i=1:length(FoldersList)
        CurrentFolder = FoldersList{i};
        ModifiedFolder = regexprep(CurrentFolder, 'pt', '.');
        numbers = regexp(ModifiedFolder, '\d+(\.\d+)?', 'match');
        numeric_values = str2double(numbers);
        Gamma = numeric_values(1);
        Wo = numeric_values(2);
        if Gamma==0.5
            L = Lengths(1);
        elseif Gamma==1.0
            L = Lengths(2);
        elseif Gamma==1.5
            L = Lengths(3);
        end
    
        p0 = 7584.3;
         h0 = 4.5e-5;
         b = 4.5e-5; 
         E = 5.9e+5;
         nu_s = 0.45;
         omega = 17951.94;%5714*pi;
         epsilon = h0/L;
         G=E/(2*(1+nu_s));
         K=E/(3*(1-2*nu_s));
         lamda = K-2/3*G;
         kappa_param = b/((lamda+2*G));
         beta = kappa_param*p0/h0;
         defscale = beta*h0;
        C_W = (1-2*nu_s)/(2*(1-nu_s))*b/G;
        C_I = (2*nu_s)*(nu_s-1/4)/(3*(1-nu_s)^2)*b^3/G;
        theta = C_I/(C_W*L^2);
    
        %% 
        pmean_files = dir(append(FoldersDir,CurrentFolder,'/pmean*.csv'));
        Z = linspace(0,1,500);
        AllData = [];
    
        for iii=1:350
            T = readtable(append(FoldersDir,CurrentFolder,'/',pmean_files(iii).name));
            AllData = [AllData,T{:,:}];  
        end
        TempAvgP = mean(AllData,2);
        fig1=figure(1);
        colorarray_sim=viridis(length(FoldersList));
        count_sim = i;
        set(fig1,'DefaultAxesFontName','Times New Roman')
        set(fig1,'DefaultAxesFontSize',14)
        h(i)=plot(Z(1:15:end),TempAvgP(1:15:end)/p0/(beta),'Square','Color',colorarray_sim(count_sim,:), 'MarkerSize',3,'MarkerFaceColor',colorarray_sim(count_sim,:));
        set(gca,'defaulttextinterpreter','latex')
        set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
        set(gca,'box','on');
        xlabel('$Z$');
        ylabel('$\langle P_1 \rangle$');
    
        hold on;
    
        %%
        syms kappa
        syms Womer Gm Th Y Z_ZR Tm real positive
        
        control = "pressure";
        Q1at0 = 0; 0.0019;
        
        % leading-order pressure
        if control == "pressure"
            P0a = sinh(kappa*(1-Z_ZR))/sinh(kappa);
        elseif control == "flow"
            P0a = sinh(kappa*(1-Z_ZR))/(kappa*cosh(kappa));
        else
            error(' ERROR: incorrect control type, flow or pressure');
        end
        P0 = P0a*exp(1i*Tm);
        
        % leading-order axial velocity
        Vz0a = 1/(1i*Womer^2)*(1 - cos(1i^(3/2)*(1-2*Y)*Womer/2)/cos(1i^(3/2)*Womer/2))*(-diff(P0a,Z_ZR));
        Vz0 = Vz0a*exp(1i*Tm);
        
        % leading-order vertical velocity, found from COM
        syms f(Y)
        ode = diff(f,Y) + diff(Vz0a,Z_ZR) == 0;
        cond = f(0) == 0;
        Vy0a = dsolve(ode,cond);
        Vy0 = Vy0a*exp(1i*Tm);
        
        % pre-compute the advective streaming terms
        inertia_avg_term1 = simplify(1/2*real(conj(Vy0)*diff(Vz0,Y)));
        inertia_avg_term2 = simplify(1/2*real(conj(Vz0)*diff(Vz0,Z_ZR)));
        
        % now define the auxilliary function -- if we defined them earlier it
        % makes the symbolic calculation harder
        f0 = (-1j*Womer - 2*1j^(3/2)*tan(1j^(3/2)*Womer/2))/Womer^3;
        kappa = sqrt(1j*Gm/(f0 + Th*1j*Gm));
        % substitute kappa and f0 into the inertia_avg_term
        inertia_avg_term1 = subs(inertia_avg_term1);
        inertia_avg_term2 = subs(inertia_avg_term2);  
        
        % numerically  
        opts = bvpset('RelTol',1e-6,'AbsTol',1e-6,'Stats','on');
        ymesh = linspace(0,1,100);
        solinit = bvpinit(ymesh, @(x)([x*(1-x); 1-2*x]));
        bcfun = @(ya,yb)([ya(1); yb(1)]);
        
        % need to set these for definiteness
        Gm = Gamma;
        Womer = Wo;
        Th = theta;
        
        % now solve for int <V_{z,1}^adv>(Y,Z) dY at each Z in zmesh
        zmesh = linspace(0,1,50);
        yint_Vz1_adv = 0*zmesh;
        for ii = 1:length(zmesh)
            % get current Z value and evaluate inertial_avg_termws
            Z_ZR = zmesh(ii);
            inertia_avg_term1_2 = subs(Womer^2/Gm*inertia_avg_term1);
            inertia_avg_term2_2 = subs(Womer^2/Gm*inertia_avg_term2);
            
            % define ODE RHS, taking into account weird behavior when 
            % argument of matlabFunction is zero
            if eval(inertia_avg_term1_2 + inertia_avg_term2_2 == 0)
                oderhs = @(Y)(0*Y);
            else
                oderhs = matlabFunction(inertia_avg_term1_2 + inertia_avg_term2_2);
            end
            bvpfun = @(Y,y)([y(2); oderhs(Y)]); 
            
            % find solution to BVP for <V_{z,1}^adv>(Y) at this Z location
            sol = bvp4c(bvpfun,bcfun,solinit,opts);
            
            % integrate from 0 to 1 dY and save it
            yint_Vz1_adv(ii) = trapz(sol.x,sol.y(1,:));
        
        end
        
        % now evaluate the streaming flow rate <Q_1> = const.
        syms Z_ZR % first clear Z value
        % this is the displacement based on the combined foundation model
        U0 = P0 - Th*diff(P0,Z_ZR,2);
        % this is the "slip" velocity, <U_0 dV_{z,0}/dY> at Y = 1
        U0dVz0dY = matlabFunction(eval(subs(1/2*real(conj(U0)*subs(diff(Vz0,Y),Y,1)))));
        temp = -1/2*U0dVz0dY(zmesh) + yint_Vz1_adv;
        
        % now evaluate the streaming pressure <P_1>(Z) 
        if control == "pressure"
            % such that <P_1>(0) = <P_1>(1) = 0
            Q1 = trapz(zmesh,temp);
            P1 = 12*cumtrapz(zmesh,temp) - 12*zmesh*Q1;
        elseif control == "flow"
            % such that <P_1>(1) = 0 and <Q_1> is given, maybe = 0
            fzmesh = flip(zmesh);
            ftemp = flip(temp);
            %Q1 = temp(end)/10*ones(size(fzmesh)); 
            Q1 = Q1at0*ones(size(fzmesh));
            %load TimeAvgQ1.mat
            %Q1 = interp1(linspace(0,1,100),norm_TempAvgQ1,fzmesh);    
            P1 = 0*ftemp;    
            % integral from 1 to Z
            for ii=1:length(zmesh)-1
                idx = length(zmesh)+1-ii;
                P1(ii) = 12*trapz(fzmesh(1:idx),ftemp(1:idx)) ...
                         + 12*(1-zmesh(ii))*Q1(idx);
            end
            P1(end) = 0;
        
        else
            error(' ERROR: incorrect control type, flow or pressure');
        end
        
        
        % % make a nice plot
        colorarray_theory=viridis(length(FoldersList));
        count_theory = i;
        set(fig1,'DefaultAxesFontName','Times New Roman')
        set(fig1,'DefaultAxesFontSize',14)
        plot(zmesh,P1,'Color',colorarray_theory(count_theory,:));
        set(gca,'defaulttextinterpreter','latex')
        set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
        set(gca,'box','on');
        xlabel('$Z$');
        ylabel('$\langle P_1 \rangle$');
        hold on;
        Gammas(i) = Gamma;
    end
    legendEntries = arrayfun(@(gm) ['$\mathrm{\gamma} = ' num2str(gm) '$'], Gammas, 'UniformOutput', false);
    legend(h,legendEntries, 'Interpreter', 'latex');
    WoStr = sprintf('%.1f', Wo);
    exportgraphics(figure(1), sprintf('/home/tools/a/urade/Downloads/Feb24/Gms_Wo_%s.pdf', WoStr), ...
         'ContentType', 'vector');   

elseif control_new=="P0"
    for i=1:length(FoldersList)
        CurrentFolder = FoldersList{i};
        ModifiedFolder = regexprep(CurrentFolder, 'pt', '.');
        numbers = regexp(ModifiedFolder, '\d+(\.\d+)?', 'match');
        numeric_values = str2double(numbers);
        Gamma = numeric_values(1);
        Wo = numeric_values(2);
        if Gamma==0.5
            L = Lengths(1);
        elseif Gamma==1.0
            L = Lengths(2);
        elseif Gamma==1.5
            L = Lengths(3);
        end
    
        p0 = 7584.3;
         h0 = 4.5e-5;
         b = 4.5e-5; 
         E = 5.9e+5;
         nu_s = 0.45;
         omega = 17951.94;%5714*pi;
         epsilon = h0/L;
         G=E/(2*(1+nu_s));
         K=E/(3*(1-2*nu_s));
         lamda = K-2/3*G;
         kappa_param = b/((lamda+2*G));
         beta = kappa_param*p0/h0;
         defscale = beta*h0;
        C_W = (1-2*nu_s)/(2*(1-nu_s))*b/G;
        C_I = (2*nu_s)*(nu_s-1/4)/(3*(1-nu_s)^2)*b^3/G;
        theta = C_I/(C_W*L^2);
    
        %% 
        pmean_files = dir(append(FoldersDir,CurrentFolder,'/pmean*.csv'));
        Z = linspace(0,1,500);
        AllData = [];
        cycle_init_index = 1;
        ref_indices = [0 70 140 210 280];
        indices_needed = cycle_init_index + ref_indices;
        fig1=figure(1);
        set(fig1,'DefaultAxesFontName','Times New Roman')
        set(fig1,'DefaultAxesFontSize',18)
        colorarray_sim=hsv(length(indices_needed));
        for iii=1:length(ref_indices)
            T = readtable(append(FoldersDir,CurrentFolder,'/',pmean_files(indices_needed(iii)).name));
            Data = T{:,:};
            count_sim = iii;
            h(iii)=plot(Z(1:15:end),Data(1:15:end)/(p0),'Square','Color',colorarray_sim(count_sim,:), 'MarkerSize',7,'MarkerFaceColor',colorarray_sim(count_sim,:));
            hold on
        end

        f0 = (1/(1i*Wo^2))*(1- (tan((1i)^1.5*Wo/2))/((1i)^1.5*Wo/2));
        k_Wo = sqrt(1i*Gamma/(f0+theta*1i*Gamma));
        omegat = omega*1e-6*ref_indices;
        colorarray_theory=hsv(length(indices_needed));  
        for iii=1:length(omegat)
            P_0 = real((sinh(k_Wo*(1-Z))/sinh(k_Wo))*exp(1i*omegat(iii)));
            count_theory = iii;
             plot(Z,P_0,'Color',colorarray_theory(count_theory,:),'LineWidth',1.5);
             hold on;
        end
        set(gca,'defaulttextinterpreter','latex')
        set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
        set(gca,'box','on');
        xlabel('$Z$');
        ylabel('$P$');
        legendEntries = {'$T =0$', '$T = 2\pi/5$', '$T = 4\pi/5$', ...
                 '$T = 6\pi/5$', '$T = 8\pi/5$'};
        leg=legend(h,legendEntries, 'Interpreter', 'latex');
        set(leg, 'FontSize', 10);
        WoStr = sprintf('%.1f', Wo);
        GmStr = sprintf('%.1f', Gamma);
        exportgraphics(figure(1), sprintf('/home/tools/a/urade/Downloads/Feb24/Gm_%s_Wo_%s_Phases.pdf', GmStr, WoStr), ...
             'ContentType', 'vector'); 
        clf;
    end
elseif control_new == "UxUy"
    for i=1:length(FoldersList)
        CurrentFolder = FoldersList{i};
        ModifiedFolder = regexprep(CurrentFolder, 'pt', '.');
        numbers = regexp(ModifiedFolder, '\d+(\.\d+)?', 'match');
        numeric_values = str2double(numbers);
        Gamma = numeric_values(1);
        Wo = numeric_values(2);
        if Gamma==0.5
            L = Lengths(1);
        elseif Gamma==1.0
            L = Lengths(2);
        elseif Gamma==1.5
            L = Lengths(3);
        end
    
        p0 = 7584.3;
         h0 = 4.5e-5;
         b = 4.5e-5; 
         E = 5.9e+5;
         nu_s = 0.45;
         omega = 17951.94;%5714*pi;
         epsilon = h0/L;
         G=E/(2*(1+nu_s));
         K=E/(3*(1-2*nu_s));
         lamda = K-2/3*G;
         kappa_param = b/((lamda+2*G));
         beta = kappa_param*p0/h0;
         defscale = beta*h0;
        C_W = (1-2*nu_s)/(2*(1-nu_s))*b/G;
        C_I = (2*nu_s)*(nu_s-1/4)/(3*(1-nu_s)^2)*b^3/G;
        theta = C_I/(C_W*L^2);
    
        %% 
        disp_files = dir(append(FoldersDir,CurrentFolder,'/disp*.csv'));
        Z = linspace(0,1,500);
        AllData = [];    
        T = readtable(append(FoldersDir,CurrentFolder,'/',disp_files((1)).name));
        AllData = [AllData,T{:,:}]; 
        disp_x = AllData(:,1);
        disp_y = AllData(:,2);
        Ux_sim = disp_x/(beta*h0);
        Uy_sim = disp_y/(beta*h0);
        fig1=figure(1);
        set(fig1,'DefaultAxesFontName','Times New Roman')
        set(fig1,'DefaultAxesFontSize',18)
        colorarray_sim=cool(2);
        count_sim=1;
        h(1)=plot(Z(1:15:end),Uy_sim(1:15:end),'o','Color','k', 'MarkerSize',7,'MarkerFaceColor','k');
        hold on;
        h(2)=plot(Z(1:15:end),Ux_sim(1:15:end),'o','Color', [0.9290 0.6940 0.1250], 'MarkerSize',7,'MarkerFaceColor',[0.9290 0.6940 0.1250]);
        hold on;

        f0 = (1/(1i*Wo^2))*(1- (tan((1i)^1.5*Wo/2))/((1i)^1.5*Wo/2));
        k_Wo = sqrt(1i*Gamma/(f0+theta*1i*Gamma));
        Uy_theory = real((f0/(1i*Gamma))*k_Wo^2*(sinh(k_Wo*(1-Z))/sinh(k_Wo))*exp(1i*0));
        colorarray_theory=cool(2);
        count_theory=1;
        plot(Z,Uy_theory,'Color','k','LineWidth',1.5);
        ylim([-0.35,1.2]);
        set(gca,'defaulttextinterpreter','latex')
        set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
        set(gca,'box','on');
        xlabel('$Z$');
        legendEntries = {'$U_{Y}$','$U_{Z}$'};
        leg=legend(h,legendEntries, 'Interpreter', 'latex');
        WoStr = sprintf('%.1f', Wo);
        GmStr = sprintf('%.1f', Gamma);
        exportgraphics(figure(1), sprintf('/home/tools/a/urade/Downloads/Feb24/Gm_%s_Wo_%s_Disp.pdf', GmStr, WoStr), ...
             'ContentType', 'vector'); 
        clf;  

    end

elseif control_new == "AvgUx"
    for i=1:length(FoldersList)
        CurrentFolder = FoldersList{i};
        ModifiedFolder = regexprep(CurrentFolder, 'pt', '.');
        numbers = regexp(ModifiedFolder, '\d+(\.\d+)?', 'match');
        numeric_values = str2double(numbers);
        Gamma = numeric_values(1);
        Wo = numeric_values(2);
        if Gamma==0.5
            L = Lengths(1);
        elseif Gamma==1.0
            L = Lengths(2);
        elseif Gamma==1.5
            L = Lengths(3);
        end
    
        p0 = 7584.3;
         h0 = 4.5e-5;
         b = 4.5e-5; 
         E = 5.9e+5;
         nu_s = 0.45;
         omega = 17951.94;%5714*pi;
         epsilon = h0/L;
         G=E/(2*(1+nu_s));
         K=E/(3*(1-2*nu_s));
         lamda = K-2/3*G;
         kappa_param = b/((lamda+2*G));
         beta = kappa_param*p0/h0;
         defscale = beta*h0;
        C_W = (1-2*nu_s)/(2*(1-nu_s))*b/G;
        C_I = (2*nu_s)*(nu_s-1/4)/(3*(1-nu_s)^2)*b^3/G;
        theta = C_I/(C_W*L^2);
    
        %% 
        disp_files = dir(append(FoldersDir,CurrentFolder,'/disp*.csv'));
        Z = linspace(0,1,500);
        AllData = [];    
        for iii=1:350
            T = readtable(append(FoldersDir,CurrentFolder,'/',disp_files(iii).name));
            AllData = [AllData,T{:,1}];  
        end
        TempAvgDispX = mean(AllData,2)/(beta^2*h0);
        Wos(i)=Wo;
        fig1=figure(1);
        set(fig1,'DefaultAxesFontName','Times New Roman')
        set(fig1,'DefaultAxesFontSize',18)
        colorarray_sim=viridis(length(FoldersList));
        count_sim=i;
        plot(Z(1:15:end),TempAvgDispX(1:15:end),'o','Color',colorarray_sim(count_sim,:), 'MarkerSize',7,'MarkerFaceColor',colorarray_sim(count_sim,:));
        hold on;

    end
        set(gca,'defaulttextinterpreter','latex')
        set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
        set(gca,'box','on');
        xlabel('$Z$');
        ylabel('$\langle U_{Z} \rangle / \beta$');
        ylim([-0.5,2.5])
        legendEntries = arrayfun(@(wo) ['$\mathrm{Wo} = ' num2str(wo) '$'], Wos, 'UniformOutput', false);
        leg=legend(legendEntries, 'Interpreter', 'latex');
        set(leg, 'Location', 'south');
        leg.Position(1) = 0.5 - 0.25; 
        leg.Position(2) = 0.22; 
        set(leg, 'FontSize', 14);
        GmStr = sprintf('%.1f', Gamma);
        exportgraphics(figure(1), sprintf('/home/tools/a/urade/Downloads/Feb24/Gm_%s_Wos_AvgDisp.pdf', GmStr), ...
             'ContentType', 'vector'); 

elseif control_new=="Pressure"
    % Initialize global min and max
    cmin = Inf;
    cmax = -Inf;
    
    % First Pass: Determine cmin and cmax across all datasets
    for i = 1:length(FoldersList)
        CurrentFolder = FoldersList{i};
            ModifiedFolder = regexprep(CurrentFolder, 'pt', '.');
        numbers = regexp(ModifiedFolder, '\d+(\.\d+)?', 'match');
        numeric_values = str2double(numbers);
        Gamma = numeric_values(1);
        Wo = numeric_values(2);
        if Gamma==0.5
            L = Lengths(1);
        elseif Gamma==1.0
            L = Lengths(2);
        elseif Gamma==1.5
            L = Lengths(3);
        end
    
        p0 = 7584.3;
         h0 = 4.5e-5;
         b = 4.5e-5; 
         E = 5.9e+5;
         nu_s = 0.45;
         omega = 17951.94;%5714*pi;
         epsilon = h0/L;
         G=E/(2*(1+nu_s));
         K=E/(3*(1-2*nu_s));
         lamda = K-2/3*G;
         kappa_param = b/((lamda+2*G));
         beta = kappa_param*p0/h0;
         defscale = beta*h0;
        C_W = (1-2*nu_s)/(2*(1-nu_s))*b/G;
        C_I = (2*nu_s)*(nu_s-1/4)/(3*(1-nu_s)^2)*b^3/G;
        theta = C_I/(C_W*L^2);    
        % Extract Pressure Data
        p_file = dir(append(FoldersDir, CurrentFolder, '/Pressure.csv'));
        data = readmatrix(append(FoldersDir, CurrentFolder, '/', p_file.name));
        
        % Extract pressure column
        pressure = data(:,1) / p0;  % Normalize by p0
        
        % Update global min and max
        cmin = min(cmin, min(pressure));
        cmax = max(cmax, max(pressure));
    end

    for i=1:length(FoldersList)
        CurrentFolder = FoldersList{i};
        ModifiedFolder = regexprep(CurrentFolder, 'pt', '.');
        numbers = regexp(ModifiedFolder, '\d+(\.\d+)?', 'match');
        numeric_values = str2double(numbers);
        Gamma = numeric_values(1);
        Wo = numeric_values(2);
        if Gamma==0.5
            L = Lengths(1);
        elseif Gamma==1.0
            L = Lengths(2);
        elseif Gamma==1.5
            L = Lengths(3);
        end
    
        p0 = 7584.3;
         h0 = 4.5e-5;
         b = 4.5e-5; 
         E = 5.9e+5;
         nu_s = 0.45;
         omega = 17951.94;%5714*pi;
         epsilon = h0/L;
         G=E/(2*(1+nu_s));
         K=E/(3*(1-2*nu_s));
         lamda = K-2/3*G;
         kappa_param = b/((lamda+2*G));
         beta = kappa_param*p0/h0;
         defscale = beta*h0;
        C_W = (1-2*nu_s)/(2*(1-nu_s))*b/G;
        C_I = (2*nu_s)*(nu_s-1/4)/(3*(1-nu_s)^2)*b^3/G;
        theta = C_I/(C_W*L^2);
    
        %% 
        p_file = dir(append(FoldersDir,CurrentFolder,'/Pressure.csv'));
        fig1=figure(1);
        set(fig1,'DefaultAxesFontName','Times New Roman')
        set(fig1,'DefaultAxesFontSize',18)

        data = readmatrix(append(FoldersDir,CurrentFolder,'/',p_file.name));
        
        % Extract columns (assuming format: X, Y, Pressure)
        x = data(:,2); 
        y = data(:,3);
        pressure = data(:,1);  % Assuming pressure is in the last column
        
        x = x/L;
        y = y/h0;
        pressure = pressure/p0;
        %% Create a Grid for Interpolation
        gridSize = 500;  % Number of interpolation points in each direction
        [Xq, Yq] = meshgrid(linspace(min(x), max(x), gridSize), linspace(min(y), max(y), gridSize));
        Pq = griddata(x, y, pressure, Xq, Yq, 'cubic'); % Interpolation using cubic method
        
        %% Plot 1: Image-style Color Plot (Similar to ParaView)
        imagesc([min(x) max(x)], [min(y) max(y)], Pq);
        axis xy; % Ensure correct axis orientation
        c = colorbar;
        colormap(parula); % Use ParaView-like color scheme
        caxis([cmin cmax]);
        set(gca,'defaulttextinterpreter','latex')
        set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
        set(gca,'box','on');
        xlabel('$Z$');
        ylabel('$Y$');
        ylabel(c,'$P$','Interpreter','latex','FontSize',18,'FontName','Times New Roman');
        WoStr = sprintf('%.1f', Wo);
        GmStr = sprintf('%.1f', Gamma);
        exportgraphics(figure(1), sprintf('/home/tools/a/urade/Downloads/Feb24/Gm_%s_Wo_%s_Pressure.pdf', GmStr, WoStr), ...
             'ContentType', 'vector'); 
        clf;
    end

elseif control_new=="Velocity"
    % Initialize global min and max
    cmin = Inf;
    cmax = -Inf;
    for i=1:length(FoldersList)
        CurrentFolder = FoldersList{i};
        ModifiedFolder = regexprep(CurrentFolder, 'pt', '.');
        numbers = regexp(ModifiedFolder, '\d+(\.\d+)?', 'match');
        numeric_values = str2double(numbers);
        Gamma = numeric_values(1);
        Wo = numeric_values(2);
        if Gamma==0.5
            L = Lengths(1);
        elseif Gamma==1.0
            L = Lengths(2);
        elseif Gamma==1.5
            L = Lengths(3);
        end
    
        p0 = 7584.3;
         h0 = 4.5e-5;
         b = 4.5e-5; 
         E = 5.9e+5;
         nu_s = 0.45;
         mu_f = 0.12;
         omega = 17951.94;%5714*pi;
         epsilon = h0/L;
         G=E/(2*(1+nu_s));
         K=E/(3*(1-2*nu_s));
         lamda = K-2/3*G;
         kappa_param = b/((lamda+2*G));
         beta = kappa_param*p0/h0;
         defscale = beta*h0;
        C_W = (1-2*nu_s)/(2*(1-nu_s))*b/G;
        C_I = (2*nu_s)*(nu_s-1/4)/(3*(1-nu_s)^2)*b^3/G;
        theta = C_I/(C_W*L^2);
        v_c = epsilon*h0*p0/mu_f;
    
        %% 
        v_file = dir(append(FoldersDir,CurrentFolder,'/Velocity.csv'));

        data = readmatrix(append(FoldersDir,CurrentFolder,'/',v_file.name));
        
        velocity_x = data(:,1);  % Assuming pressure is in the last column
        velocity_x = velocity_x/v_c;
        % Update global min and max
        cmin = min(cmin, min(velocity_x));
        cmax = max(cmax, max(velocity_x));
    end

    for i=1:length(FoldersList)
        CurrentFolder = FoldersList{i};
        ModifiedFolder = regexprep(CurrentFolder, 'pt', '.');
        numbers = regexp(ModifiedFolder, '\d+(\.\d+)?', 'match');
        numeric_values = str2double(numbers);
        Gamma = numeric_values(1);
        Wo = numeric_values(2);
        if Gamma==0.5
            L = Lengths(1);
        elseif Gamma==1.0
            L = Lengths(2);
        elseif Gamma==1.5
            L = Lengths(3);
        end
    
        p0 = 7584.3;
         h0 = 4.5e-5;
         b = 4.5e-5; 
         E = 5.9e+5;
         nu_s = 0.45;
         mu_f = 0.12;
         omega = 17951.94;%5714*pi;
         epsilon = h0/L;
         G=E/(2*(1+nu_s));
         K=E/(3*(1-2*nu_s));
         lamda = K-2/3*G;
         kappa_param = b/((lamda+2*G));
         beta = kappa_param*p0/h0;
         defscale = beta*h0;
        C_W = (1-2*nu_s)/(2*(1-nu_s))*b/G;
        C_I = (2*nu_s)*(nu_s-1/4)/(3*(1-nu_s)^2)*b^3/G;
        theta = C_I/(C_W*L^2);
        v_c = epsilon*h0*p0/mu_f;
    
        %% 
        v_file = dir(append(FoldersDir,CurrentFolder,'/Velocity.csv'));
        fig1=figure(1);
        set(fig1,'DefaultAxesFontName','Times New Roman')
        set(fig1,'DefaultAxesFontSize',18)

        data = readmatrix(append(FoldersDir,CurrentFolder,'/',v_file.name));
        
        % Extract columns (assuming format: X, Y, Pressure)
        x = data(:,4); 
        y = data(:,5);
        velocity_x = data(:,1);  % Assuming pressure is in the last column
        
        x = x/L;
        y = y/h0;
        velocity_x = velocity_x/v_c;
        %% Create a Grid for Interpolation
        gridSize = 500;  % Number of interpolation points in each direction
        [Xq, Yq] = meshgrid(linspace(min(x), max(x), gridSize), linspace(min(y), max(y), gridSize));
        Vq = griddata(x, y, velocity_x, Xq, Yq, 'cubic'); % Interpolation using cubic method
        
        %% Plot 1: Image-style Color Plot (Similar to ParaView)
        imagesc([min(x) max(x)], [min(y) max(y)], Vq);
        axis xy; % Ensure correct axis orientation
        c = colorbar;
        colormap(jet); % Use ParaView-like color scheme
        caxis([cmin cmax]);
        set(gca,'defaulttextinterpreter','latex')
        set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
        set(gca,'box','on');
        xlabel('$Z$');
        ylabel('$Y$');
        ylabel(c,'$V_Z$','Interpreter','latex','FontSize',18,'FontName','Times New Roman');
            tickValues = linspace(min(c.Limits), max(c.Limits), 3);  % Set 5 evenly spaced tick marks
            c.Ticks = tickValues;  
            c.TickLabels = compose('%.1f', tickValues);  % Format to one decimal place
        WoStr = sprintf('%.1f', Wo);
        GmStr = sprintf('%.1f', Gamma);
        exportgraphics(figure(1), sprintf('/home/tools/a/urade/Downloads/Feb24/Gm_%s_Wo_%s_Vx.pdf', GmStr, WoStr), ...
            'ContentType', 'vector'); 
        clf;
    end

elseif control_new=="Mesh"
    % Initialize global min and max
    cmin_x = Inf;
    cmax_x = -Inf;
    cmin_y = Inf;
    cmax_y = -Inf;
    for i=1:length(FoldersList)
        CurrentFolder = FoldersList{i};
        ModifiedFolder = regexprep(CurrentFolder, 'pt', '.');
        numbers = regexp(ModifiedFolder, '\d+(\.\d+)?', 'match');
        numeric_values = str2double(numbers);
        Gamma = numeric_values(1);
        Wo = numeric_values(2);
        if Gamma==0.5
            L = Lengths(1);
        elseif Gamma==1.0
            L = Lengths(2);
        elseif Gamma==1.5
            L = Lengths(3);
        end
    
        p0 = 7584.3;
         h0 = 4.5e-5;
         b = 4.5e-5; 
         E = 5.9e+5;
         nu_s = 0.45;
         mu_f = 0.12;
         omega = 17951.94;%5714*pi;
         epsilon = h0/L;
         G=E/(2*(1+nu_s));
         K=E/(3*(1-2*nu_s));
         lamda = K-2/3*G;
         kappa_param = b/((lamda+2*G));
         beta = kappa_param*p0/h0;
         defscale = beta*h0;
        C_W = (1-2*nu_s)/(2*(1-nu_s))*b/G;
        C_I = (2*nu_s)*(nu_s-1/4)/(3*(1-nu_s)^2)*b^3/G;
        theta = C_I/(C_W*L^2);
        v_c = epsilon*h0*p0/mu_f;
    
        %% 
        m_file = dir(append(FoldersDir,CurrentFolder,'/Mesh.csv'));

        data = readmatrix(append(FoldersDir,CurrentFolder,'/',m_file.name));
        
        % Extract columns (assuming format: X, Y, Pressure)
        x = data(:,4); 
        y = data(:,5);
        mesh_x = data(:,1);  % Assuming pressure is in the last column
        mesh_y = data(:,2); 
        x = x/L;
        y = y/h0;
        mesh_x = mesh_x/(beta*h0);
        mesh_y = mesh_y/(beta*h0);
        idx = (y>=1);
        x = x(idx);
        y = y(idx);
        mesh_x = mesh_x(idx);
        mesh_y = mesh_y(idx);
        % Update global min and max
        cmin_x = min(cmin_x, min(mesh_x));
        cmax_x = max(cmax_x, max(mesh_x));
        cmin_y = min(cmin_y, min(mesh_y));
        cmax_y = max(cmax_y, max(mesh_y));
    end
    for i=1:length(FoldersList)
        CurrentFolder = FoldersList{i};
        ModifiedFolder = regexprep(CurrentFolder, 'pt', '.');
        numbers = regexp(ModifiedFolder, '\d+(\.\d+)?', 'match');
        numeric_values = str2double(numbers);
        Gamma = numeric_values(1);
        Wo = numeric_values(2);
        if Gamma==0.5
            L = Lengths(1);
        elseif Gamma==1.0
            L = Lengths(2);
        elseif Gamma==1.5
            L = Lengths(3);
        end
    
        p0 = 7584.3;
         h0 = 4.5e-5;
         b = 4.5e-5; 
         E = 5.9e+5;
         nu_s = 0.45;
         mu_f = 0.12;
         omega = 17951.94;%5714*pi;
         epsilon = h0/L;
         G=E/(2*(1+nu_s));
         K=E/(3*(1-2*nu_s));
         lamda = K-2/3*G;
         kappa_param = b/((lamda+2*G));
         beta = kappa_param*p0/h0;
         defscale = beta*h0;
        C_W = (1-2*nu_s)/(2*(1-nu_s))*b/G;
        C_I = (2*nu_s)*(nu_s-1/4)/(3*(1-nu_s)^2)*b^3/G;
        theta = C_I/(C_W*L^2);
        v_c = epsilon*h0*p0/mu_f;
    
        %% 
        m_file = dir(append(FoldersDir,CurrentFolder,'/Mesh.csv'));
        fig1=figure(1);
        set(fig1,'DefaultAxesFontName','Times New Roman')
        set(fig1,'DefaultAxesFontSize',18)

        data = readmatrix(append(FoldersDir,CurrentFolder,'/',m_file.name));
        
        % Extract columns (assuming format: X, Y, Pressure)
        x = data(:,4); 
        y = data(:,5);
        mesh_x = data(:,1);  % Assuming pressure is in the last column
        mesh_y = data(:,2); 
        x = x/L;
        y = y/h0;
        mesh_x = mesh_x/(beta*h0);
        mesh_y = mesh_y/(beta*h0);
        idx = (y>=1);
        x = x(idx);
        y = y(idx);
        mesh_x = mesh_x(idx);
        mesh_y = mesh_y(idx);       
        %% Create a Grid for Interpolation
        gridSize = 500;  % Number of interpolation points in each direction
        [Xq, Yq] = meshgrid(linspace(min(x), max(x), gridSize), linspace(min(y), max(y), gridSize));
        Mqx = griddata(x, y, mesh_x, Xq, Yq, 'cubic'); % Interpolation using cubic method
        Mqy = griddata(x, y, mesh_y, Xq, Yq, 'cubic'); % Interpolation using cubic method
        %% Plot 1: Image-style Color Plot (Similar to ParaView)
        imagesc([min(x) max(x)], [min(y) max(y)], Mqx);
        axis xy; % Ensure correct axis orientation
        c = colorbar;
        colormap(jet); % Use ParaView-like color scheme
        caxis([cmin_x cmax_x]);
        set(gca,'defaulttextinterpreter','latex')
        set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
        set(gca,'box','on');
        xlabel('$Z$');
        ylabel('$Y$');
        ylabel(c,'$U_Z$','Interpreter','latex','FontSize',18,'FontName','Times New Roman');
            tickValues = linspace(min(c.Limits), max(c.Limits), 3);  % Set 5 evenly spaced tick marks
            c.Ticks = tickValues;  
            c.TickLabels = compose('%.1f', tickValues);  % Format to one decimal place
        WoStr = sprintf('%.1f', Wo);
        GmStr = sprintf('%.1f', Gamma);
        exportgraphics(figure(1), sprintf('/home/tools/a/urade/Downloads/Feb24/Gm_%s_Wo_%s_Ux.pdf', GmStr, WoStr), ...
            'ContentType', 'vector'); 
        clf;

        %% Plot 2: Image-style Color Plot (Similar to ParaView)
        imagesc([min(x) max(x)], [min(y) max(y)], Mqy);
        axis xy; % Ensure correct axis orientation
        c = colorbar;
        colormap(jet); % Use ParaView-like color scheme
        caxis([cmin_y cmax_y]);
        set(gca,'defaulttextinterpreter','latex')
        set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
        set(gca,'box','on');
        xlabel('$Z$');
        ylabel('$Y$');
        ylabel(c,'$U_Y$','Interpreter','latex','FontSize',18,'FontName','Times New Roman');
            tickValues = linspace(min(c.Limits), max(c.Limits), 3);  % Set 5 evenly spaced tick marks
            c.Ticks = tickValues;  
            c.TickLabels = compose('%.1f', tickValues);  % Format to one decimal place
        WoStr = sprintf('%.1f', Wo);
        GmStr = sprintf('%.1f', Gamma);
        exportgraphics(figure(1), sprintf('/home/tools/a/urade/Downloads/Feb24/Gm_%s_Wo_%s_Uy.pdf', GmStr, WoStr), ...
            'ContentType', 'vector'); 
        clf;
    end

elseif control_new=="AvgUy"
    for i=1:length(FoldersList)
        CurrentFolder = FoldersList{i};
        ModifiedFolder = regexprep(CurrentFolder, 'pt', '.');
        numbers = regexp(ModifiedFolder, '\d+(\.\d+)?', 'match');
        numeric_values = str2double(numbers);
        Gamma = numeric_values(1);
        Wo = numeric_values(2);
        if Gamma==0.5
            L = Lengths(1);
        elseif Gamma==1.0
            L = Lengths(2);
        elseif Gamma==1.5
            L = Lengths(3);
        end
    
        p0 = 7584.3;
         h0 = 4.5e-5;
         b = 4.5e-5; 
         E = 5.9e+5;
         nu_s = 0.45;
         omega = 17951.94;%5714*pi;
         epsilon = h0/L;
         G=E/(2*(1+nu_s));
         K=E/(3*(1-2*nu_s));
         lamda = K-2/3*G;
         kappa_param = b/((lamda+2*G));
         beta = kappa_param*p0/h0;
         defscale = beta*h0;
        C_W = (1-2*nu_s)/(2*(1-nu_s))*b/G;
        C_I = (2*nu_s)*(nu_s-1/4)/(3*(1-nu_s)^2)*b^3/G;
        theta = C_I/(C_W*L^2);
    
        %% 
        disp_files = dir(append(FoldersDir,CurrentFolder,'/disp*.csv'));
        Z = linspace(0,1,500);
        AllData = [];
        %time = linspace(0,2*pi,351);
        for iii=1:350
            T = readtable(append(FoldersDir,CurrentFolder,'/',disp_files(iii).name));
            AllData = [AllData,T{:,2}];  
        end
        
        % for iii=1:500
        %     TempAvgUy(iii)=trapz(time,AllData(iii,:))/(time(end)-time(1))/(beta^2*h0);
        % end
        TempAvgUy = mean(AllData,2)/(beta^2*h0);
        fig1=figure(1);
        colorarray_sim=viridis(length(FoldersList));
        count_sim = i;
        set(fig1,'DefaultAxesFontName','Times New Roman')
        set(fig1,'DefaultAxesFontSize',18)
        h(i)=plot(Z(1:15:end),TempAvgUy(1:15:end),'o','Color',colorarray_sim(count_sim,:), 'MarkerSize',7,'MarkerFaceColor',colorarray_sim(count_sim,:));
        set(gca,'defaulttextinterpreter','latex')
        set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
        set(gca,'box','on');
        xlabel('$Z$');
        ylabel('$\langle U_Y \rangle / \beta$');
    
        hold on;
    
        %%
        syms kappa
        syms Womer Gm Th Y Z_ZR Tm real positive
        
        control = "pressure";
        Q1at0 = 0; 0.0019;
        
        % leading-order pressure
        if control == "pressure"
            P0a = sinh(kappa*(1-Z_ZR))/sinh(kappa);
        elseif control == "flow"
            P0a = sinh(kappa*(1-Z_ZR))/(kappa*cosh(kappa));
        else
            error(' ERROR: incorrect control type, flow or pressure');
        end
        P0 = P0a*exp(1i*Tm);
        
        % leading-order axial velocity
        Vz0a = 1/(1i*Womer^2)*(1 - cos(1i^(3/2)*(1-2*Y)*Womer/2)/cos(1i^(3/2)*Womer/2))*(-diff(P0a,Z_ZR));
        Vz0 = Vz0a*exp(1i*Tm);
        
        % leading-order vertical velocity, found from COM
        syms f(Y)
        ode = diff(f,Y) + diff(Vz0a,Z_ZR) == 0;
        cond = f(0) == 0;
        Vy0a = dsolve(ode,cond);
        Vy0 = Vy0a*exp(1i*Tm);
        
        % pre-compute the advective streaming terms
        inertia_avg_term1 = simplify(1/2*real(conj(Vy0)*diff(Vz0,Y)));
        inertia_avg_term2 = simplify(1/2*real(conj(Vz0)*diff(Vz0,Z_ZR)));
        
        % now define the auxilliary function -- if we defined them earlier it
        % makes the symbolic calculation harder
        f0 = (-1j*Womer - 2*1j^(3/2)*tan(1j^(3/2)*Womer/2))/Womer^3;
        kappa = sqrt(1j*Gm/(f0 + Th*1j*Gm));
        % substitute kappa and f0 into the inertia_avg_term
        inertia_avg_term1 = subs(inertia_avg_term1);
        inertia_avg_term2 = subs(inertia_avg_term2);  
        
        % numerically  
        opts = bvpset('RelTol',1e-6,'AbsTol',1e-6,'Stats','on');
        ymesh = linspace(0,1,100);
        solinit = bvpinit(ymesh, @(x)([x*(1-x); 1-2*x]));
        bcfun = @(ya,yb)([ya(1); yb(1)]);
        
        % need to set these for definiteness
        Gm = Gamma;
        Womer = Wo;
        Th = theta;
        
        % now solve for int <V_{z,1}^adv>(Y,Z) dY at each Z in zmesh
        zmesh = linspace(0,1,50);
        yint_Vz1_adv = 0*zmesh;
        for ii = 1:length(zmesh)
            % get current Z value and evaluate inertial_avg_termws
            Z_ZR = zmesh(ii);
            inertia_avg_term1_2 = subs(Womer^2/Gm*inertia_avg_term1);
            inertia_avg_term2_2 = subs(Womer^2/Gm*inertia_avg_term2);
            
            % define ODE RHS, taking into account weird behavior when 
            % argument of matlabFunction is zero
            if eval(inertia_avg_term1_2 + inertia_avg_term2_2 == 0)
                oderhs = @(Y)(0*Y);
            else
                oderhs = matlabFunction(inertia_avg_term1_2 + inertia_avg_term2_2);
            end
            bvpfun = @(Y,y)([y(2); oderhs(Y)]); 
            
            % find solution to BVP for <V_{z,1}^adv>(Y) at this Z location
            sol = bvp4c(bvpfun,bcfun,solinit,opts);
            
            % integrate from 0 to 1 dY and save it
            yint_Vz1_adv(ii) = trapz(sol.x,sol.y(1,:));
        
        end
        
        % now evaluate the streaming flow rate <Q_1> = const.
        syms Z_ZR % first clear Z value
        % this is the displacement based on the combined foundation model
        U0 = P0 - Th*diff(P0,Z_ZR,2);
        % this is the "slip" velocity, <U_0 dV_{z,0}/dY> at Y = 1
        U0dVz0dY = matlabFunction(eval(subs(1/2*real(conj(U0)*subs(diff(Vz0,Y),Y,1)))));
        temp = -1/2*U0dVz0dY(zmesh) + yint_Vz1_adv;
        
        % now evaluate the streaming pressure <P_1>(Z) 
        if control == "pressure"
            % such that <P_1>(0) = <P_1>(1) = 0
            Q1 = trapz(zmesh,temp);
            P1 = 12*cumtrapz(zmesh,temp) - 12*zmesh*Q1;
        elseif control == "flow"
            % such that <P_1>(1) = 0 and <Q_1> is given, maybe = 0
            fzmesh = flip(zmesh);
            ftemp = flip(temp);
            %Q1 = temp(end)/10*ones(size(fzmesh)); 
            Q1 = Q1at0*ones(size(fzmesh));
            %load TimeAvgQ1.mat
            %Q1 = interp1(linspace(0,1,100),norm_TempAvgQ1,fzmesh);    
            P1 = 0*ftemp;    
            % integral from 1 to Z
            for ii=1:length(zmesh)-1
                idx = length(zmesh)+1-ii;
                P1(ii) = 12*trapz(fzmesh(1:idx),ftemp(1:idx)) ...
                         + 12*(1-zmesh(ii))*Q1(idx);
            end
            P1(end) = 0;
        
        else
            error(' ERROR: incorrect control type, flow or pressure');
        end
        %TimeAvgUyTheory = P1-theta*gradient(gradient(P1,zmesh),zmesh);
        TimeAvgUyTheory = P1-theta*4*del2(P1)/(zmesh(2)-zmesh(1))^2;
        % % Define the Z-axis for Data1 and Data2
        % Z1 = linspace(0, 1, 50);   % Z-axis for Data1 (50 points)
        % Z2 = linspace(0, 1, 500);  % Z-axis for Data2 (500 points)
        % 
        % % Interpolate Data1 to match the Z-axis of Data2
        % Data1_interp = interp1(Z1, TimeAvgUyTheory, Z2, 'linear', 'extrap');
        % 
        % XX(i,:)=Data1_interp./TempAvgUy';

        % % Compute the spacing between mesh points (assuming uniform spacing)
        % dz = zmesh(2) - zmesh(1);
        % 
        % % Initialize the second derivative array
        % d2P1_dz2 = zeros(size(P1));
        % 
        % % Compute the second derivative using central difference
        % for j = 2:length(zmesh)-1
        %     d2P1_dz2(j) = (P1(j+1) - 2*P1(j) + P1(j-1)) / (dz^2);
        % end
        % d2P1_dz2(1) = (P1(2) - 2*P1(1)) / (dz^2);
        % d2P1_dz2(end) = (-2*P1(end)+P1(end-1)) / (dz^2);
        % TimeAvgUyTheory = P1 - theta * d2P1_dz2;
        % % make a nice plot
        colorarray_theory=viridis(length(FoldersList));
        count_theory = i;
        set(fig1,'DefaultAxesFontName','Times New Roman')
        set(fig1,'DefaultAxesFontSize',18)
        plot(zmesh,TimeAvgUyTheory,'Color',colorarray_theory(count_theory,:),'LineWidth',1.5);
        set(gca,'defaulttextinterpreter','latex')
        set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
        set(gca,'box','on');
        xlabel('$Z$');
        ylabel('$\langle U_Y \rangle / \beta$');
        hold on;
        Wos(i) = Wo;
    end
    ylim([-0.5,2.5]);
    legendEntries = arrayfun(@(wo) ['$\mathrm{Wo} = ' num2str(wo) '$'], Wos, 'UniformOutput', false);
    leg=legend(h,legendEntries, 'Interpreter', 'latex');
    % set(leg, 'Location', 'south');
    % leg.Position(1) = 0.5 - leg.Position(3)/2; 
    % leg.Position(2) = 0.22; 
    GammaStr = strrep(sprintf('%.1f', Gamma), '.', 'pt');
    exportgraphics(figure(1), sprintf('/home/tools/a/urade/Downloads/Feb24/Gm_%s_Wos_AvgUy.pdf', GammaStr), ...
         'ContentType', 'vector');  

elseif control_new=="U0"
    for i=1:length(FoldersList)
        CurrentFolder = FoldersList{i};
        ModifiedFolder = regexprep(CurrentFolder, 'pt', '.');
        numbers = regexp(ModifiedFolder, '\d+(\.\d+)?', 'match');
        numeric_values = str2double(numbers);
        Gamma = numeric_values(1);
        Wo = numeric_values(2);
        if Gamma==0.5
            L = Lengths(1);
        elseif Gamma==1.0
            L = Lengths(2);
        elseif Gamma==1.5
            L = Lengths(3);
        end
    
        p0 = 7584.3;
         h0 = 4.5e-5;
         b = 4.5e-5; 
         E = 5.9e+5;
         nu_s = 0.45;
         omega = 17951.94;%5714*pi;
         epsilon = h0/L;
         G=E/(2*(1+nu_s));
         K=E/(3*(1-2*nu_s));
         lamda = K-2/3*G;
         kappa_param = b/((lamda+2*G));
         beta = kappa_param*p0/h0;
         defscale = beta*h0;
        C_W = (1-2*nu_s)/(2*(1-nu_s))*b/G;
        C_I = (2*nu_s)*(nu_s-1/4)/(3*(1-nu_s)^2)*b^3/G;
        theta = C_I/(C_W*L^2);
    
        %% 
        disp_files = dir(append(FoldersDir,CurrentFolder,'/disp*.csv'));
        Z = linspace(0,1,500);
        AllData = [];
        cycle_init_index = 1;
        ref_indices = [0 70 140 210 280];
        indices_needed = cycle_init_index + ref_indices;
        fig1=figure(1);
        set(fig1,'DefaultAxesFontName','Times New Roman')
        set(fig1,'DefaultAxesFontSize',18)
        colorarray_sim=hsv(length(indices_needed));
        for iii=1:length(ref_indices)
            T = readtable(append(FoldersDir,CurrentFolder,'/',disp_files(indices_needed(iii)).name));
            Data = T{:,2};
            count_sim = iii;
            h(iii)=plot(Z(1:15:end),Data(1:15:end)/(beta*h0),'o','Color',colorarray_sim(count_sim,:), 'MarkerSize',7,'MarkerFaceColor',colorarray_sim(count_sim,:));
            hold on
        end

        f0 = (1/(1i*Wo^2))*(1- (tan((1i)^1.5*Wo/2))/((1i)^1.5*Wo/2));
        k_Wo = sqrt(1i*Gamma/(f0+theta*1i*Gamma));
        omegat = omega*1e-6*ref_indices;
        colorarray_theory=hsv(length(indices_needed));  
        for iii=1:length(omegat)
            U_0 = real((f0/(1i*Gamma))*k_Wo^2*(sinh(k_Wo*(1-Z))/sinh(k_Wo))*exp(1i*omegat(iii)));
            count_theory = iii;
             plot(Z,U_0,'Color',colorarray_theory(count_theory,:),'LineWidth',1.5);
             hold on;
        end
        set(gca,'defaulttextinterpreter','latex')
        set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
        set(gca,'box','on');
        xlabel('$Z$');
        ylabel('$U_Y$');
        legendEntries = {'$T =0$', '$T = 2\pi/5$', '$T = 4\pi/5$', ...
                 '$T = 6\pi/5$', '$T = 8\pi/5$'};
        leg=legend(h,legendEntries, 'Interpreter', 'latex');
        set(leg, 'FontSize', 10);
        WoStr = sprintf('%.1f', Wo);
        GmStr = sprintf('%.1f', Gamma);
        exportgraphics(figure(1), sprintf('/home/tools/a/urade/Downloads/Feb24/Gm_%s_Wo_%s_U0Phases.pdf', GmStr, WoStr), ...
             'ContentType', 'vector'); 
        clf;
    end
elseif control_new=="U0aP0aAppendix"
    for i=1:length(FoldersList)
        CurrentFolder = FoldersList{i};
        ModifiedFolder = regexprep(CurrentFolder, 'pt', '.');
        numbers = regexp(ModifiedFolder, '\d+(\.\d+)?', 'match');
        numeric_values = str2double(numbers);
        Gamma = numeric_values(1);
        Wo = numeric_values(2);
        if Gamma==0.5
            L = Lengths(1);
        elseif Gamma==1.0
            L = Lengths(2);
        elseif Gamma==1.5
            L = Lengths(3);
        end
    
        p0 = 7584.3;
         h0 = 4.5e-5;
         b = 4.5e-5; 
         E = 5.9e+5;
         nu_s = 0.45;
         omega = 17951.94;%5714*pi;
         epsilon = h0/L;
         G=E/(2*(1+nu_s));
         K=E/(3*(1-2*nu_s));
         lamda = K-2/3*G;
         kappa_param = b/((lamda+2*G));
         beta = kappa_param*p0/h0;
         defscale = beta*h0;
        C_W = (1-2*nu_s)/(2*(1-nu_s))*b/G;
        C_I = (2*nu_s)*(nu_s-1/4)/(3*(1-nu_s)^2)*b^3/G;
        theta = C_I/(C_W*L^2);
    
        %% 
        Z = linspace(0,1,500);
        AllData = [];
        TenCoeff = [10^(-3) 10^(-4) 10^(-5)];

        fig1=figure(1);
        set(fig1,'DefaultAxesFontName','Times New Roman')
        set(fig1,'DefaultAxesFontSize',18)
        fig2=figure(2);
        set(fig2,'DefaultAxesFontName','Times New Roman')
        set(fig2,'DefaultAxesFontSize',18)        
        colorarray_theory=hot(2*length(TenCoeff));

        f0 = (1/(1i*Wo^2))*(1- (tan((1i)^1.5*Wo/2))/((1i)^1.5*Wo/2));
        k_Wo = sqrt(1i*Gamma/(f0+theta*1i*Gamma));
         
        for iii=1:length(TenCoeff)
            lamda1 = sqrt((1/(2*f0*TenCoeff(iii)))*(f0+1i*Gamma*theta+sqrt(f0^2-Gamma^2*theta^2+2i*Gamma*f0*(theta-2*TenCoeff(iii)))));
            lamda2 = sqrt((1/(2*f0*TenCoeff(iii)))*(f0+1i*Gamma*theta-sqrt(f0^2-Gamma^2*theta^2+2i*Gamma*f0*(theta-2*TenCoeff(iii)))));
            a11P0a = lamda2^2/(lamda2^2-lamda1^2);
            a11U0a = a11P0a;
            P0a = real(a11P0a*(sinh(lamda1*(1-Z))/sinh(lamda1)-(lamda1/lamda2)^2*sinh(lamda2*(1-Z))/sinh(lamda2)));
            U0a = real(a11U0a*lamda1^2*f0/1i/Gamma*(sinh(lamda1*(1-Z))/sinh(lamda1)-sinh(lamda2*(1-Z))/sinh(lamda2)));
            count_theory = iii;
            figure(1);
             h(iii)=plot(Z,P0a,'Color',colorarray_theory(count_theory,:),'LineWidth',1.5);
             hold on;
            figure(2);
             g(iii)=plot(Z,U0a,'Color',colorarray_theory(count_theory,:),'LineWidth',1.5);
             hold on;
        end


        figure(1);
        set(gca,'defaulttextinterpreter','latex')
        set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
        set(gca,'box','on');
        xlabel('$Z$');
        ylabel('$\mathrm{Re}[P_{0,a}]$');
        legendEntries = {'$\mathcal{T}=10^{-3}$','$\mathcal{T}=10^{-4}$','$\mathcal{T}=10^{-5}$'};
        leg=legend(h,legendEntries, 'Interpreter', 'latex');
        set(leg, 'FontSize', 10);
        WoStr = sprintf('%.1f', Wo);
        GmStr = sprintf('%.1f', Gamma);
        exportgraphics(figure(1), sprintf('/home/tools/a/urade/Downloads/Feb24/Gm_%s_Wo_%s_P0Appen.pdf', GmStr, WoStr), ...
             'ContentType', 'vector'); 
        figure(2);
        set(gca,'defaulttextinterpreter','latex')
        set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
        set(gca,'box','on');
        xlabel('$Z$');
        ylabel('$\mathrm{Re}[U_{Y,0,a}]$');
        legendEntries = {'$\mathcal{T}=10^{-3}$','$\mathcal{T}=10^{-4}$','$\mathcal{T}=10^{-5}$'};
        leg=legend(g,legendEntries, 'Interpreter', 'latex');
        set(leg, 'FontSize', 10);
        WoStr = sprintf('%.1f', Wo);
        GmStr = sprintf('%.1f', Gamma);
        exportgraphics(figure(2), sprintf('/home/tools/a/urade/Downloads/Feb24/Gm_%s_Wo_%s_U0Appen.pdf', GmStr, WoStr), ...
             'ContentType', 'vector');         
        
    end

end


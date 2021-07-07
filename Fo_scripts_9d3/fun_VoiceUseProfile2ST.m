function fun_VoiceUseProfile2ST(F0_space,SPL_space,F0_SPL_array,filename,xlabtext,ylabtext,fignum)
% Script originally writen by Albert.  I have notated it and changed things
% for my own benifit and understanding.
%
% ejh@ncvs@dcpa 20070119
% ejh@msu 20180413
% clear
% [filename, pathname] = uigetfile('*summary*.mat', 'Pick an summary file');
% load(filename)
%%
% % load F0_SPL_matrices
% % F0_matrix is the Fo information
% % SPL_matrix is the corresponding SPL information
% %

% F0_SPL_array=round(10*F0_SPL_array)/10;
F0_SPL_space = zeros(length(SPL_space),length(F0_space)); % Fo and SPL Mtrx
SPL_array = repmat(SPL_space',1,length(F0_space)); % possible array


F0_SPL_mat_sorted = sortrows(round(F0_SPL_array));
i = 1; F0_SPL_mat = [];
for n = 2:length(F0_SPL_array)
    if F0_SPL_mat_sorted(n,2) ~= F0_SPL_mat_sorted(n-1,2)
        %         keyboard
        F0_SPL_mat = vertcat(F0_SPL_mat, [F0_SPL_mat_sorted(n-1,:) n-i]);
        i = n;
    end
end
F0_SPL_mat = vertcat(F0_SPL_mat, [F0_SPL_mat_sorted(n,:) 1]); %catch last 1
[F0_SPL_mat_max,F0_SPL_mat_maxi]=max(F0_SPL_mat);

%%
% *************************************************************************
% Creates an array for use in surface plots
% *************************************************************************


% add to be a range within two parts of fo space and spl space
for n=1:length(F0_space)-1
    
    for l=1:length(F0_SPL_array)
        if F0_SPL_array(l,1)>=F0_space(n)&&F0_SPL_array(l,1)<F0_space(n+1)
            for m=1:length(SPL_space)-1
                if F0_SPL_array(l,2)>=SPL_space(m)&&F0_SPL_array(l,2)<SPL_space(m+1)
                    F0_SPL_space(m,n)=F0_SPL_space(m,n)+1;
                end
            end
        end
    end
end
%%
% n = 1;
% for m = 1:numel(SPL_array)
%     disp([SPL_array(m) F0_SPL_mat(n,2)])
%   if SPL_array(m) == F0_SPL_mat(n,2)
%     F0_SPL_array(m) = F0_SPL_mat(n,3);
%     n = n + 1;
%   end
%   if n > length(F0_SPL_mat) break; end
% end

%%
if n>1
    %%
    % *************************************************************************
    % Set plotting parameters
    % *************************************************************************
    F0_axes = [min(F0_space) max(F0_space)];
    SPL_axes = [min(SPL_space) max(SPL_space)];
    view_angle = [45 45];
    
    
    figure(fignum+1)
    h = surf(F0_space,SPL_space,F0_SPL_space,'EdgeColor','none',...
        'FaceLighting','phong','FaceColor','interp','DiffuseStrength',1,...
        'SpecularStrength',0.3);
    xlim(F0_axes)
    %   xlim(F0_axes+[0 -600])
    ylim(SPL_axes)
    xlabel(xlabtext); ylabel(ylabtext); zlabel('# Occurrences');
    title( strrep(filename, '_', '-'))
    grid on
    view(view_angle)
    
    x = linspace(0,1,64);
    cmap = colormap(winter);
    colormap(cmap(64:-1:1,:));
    %   alpha('z');
    alpha_map = x.^0.1; % this removes most of the low count occurrences
    alphamap(alpha_map);
    lightangle(146,20)
    
%     OptionZ.FrameRate=30;OptionZ.Duration=5.5;OptionZ.Periodic=true;
%     CaptureFigVid([-30,30;-110,10;-190,80;-290,30;-400,40],...
%         [filename '_VUProtate'],OptionZ)
    
    
    %% make a contour plot
    figure(fignum+2), clf
    [C,h]=contour(F0_space,SPL_space,F0_SPL_space,15);
    xlabel(xlabtext); ylabel(ylabtext); zlabel('# Occurrences');
    hold on,
    plot(F0_SPL_mat(F0_SPL_mat_maxi(3),1),F0_SPL_mat(F0_SPL_mat_maxi(3),2),'r+'),
    % plot(F0_SPL_mat_mean(1),F0_SPL_mat_mean(2),'g+'),
    hold off
    title( strrep(filename, '_', '-'))
    
    
    %% make new VUP for contours...
    
    % get array with number of instances with the spl and fo coordinants.
    % F0_SPL_space(numinstances,Fo,SPL,Fo_space index, SPL_space index,...
    %  percent of overall in that spot, distance from maximum)
    % indx - index of F0_SPL_space sorted numinstances going from high to low
    
    totalInst=sum(sum(F0_SPL_space));
    tmp=[];
    for n=1:length(SPL_space)
        tmp=horzcat(tmp,...
            vertcat(F0_SPL_space(n,:),F0_space,F0_space.*0+SPL_space(n),...
            1:length(F0_space),F0_space.*0+n,100*F0_SPL_space(n,:)/totalInst));
        %    n=1
        %    horzcat(vertcat(F0_SPL_array(n,:),F0_space,F0_space.*0+SPL_space(n)),...
        %      vertcat(F0_SPL_array(n+1,:),F0_space,F0_space.*0+SPL_space(n+1)))
    end
    F0_SPL_space=tmp; clear tmp
    
    [B,IX] = sort(F0_SPL_space',1,'descend');
    indx=IX(1:numel(F0_SPL_space));
    clear B IX
    
    F0_SPL_space(1,indx(1));
    
    % get distance from highest or mode
    for n=1:length(F0_SPL_space)
        dist1(n)=sqrt((F0_SPL_space(2,indx(1))-F0_SPL_space(2,n))^2 + ...
            (F0_SPL_space(3,indx(1))-F0_SPL_space(3,n))^2);
    end
    
    F0_SPL_space=vertcat(F0_SPL_space,dist1); clear dist1
    
    % get cummulative sum of the percentage in order of sorted #instances
    cumsum1=100*cumsum(F0_SPL_space(1,indx))/totalInst;
    
    % find area of a single box and also get
    diffSPL=mean(diff(SPL_space))/2;
    diffF0=mean(diff(F0_space))/2;
    
    % find percentage of interest areas
    PercentsOfInterest=[68.2 50 25 10];
%     figure(fignum+3), clf, hold on
%     figure(fignum+4), clf, hold on
    %% loop through and find areas of interest
    for nn=1:length(PercentsOfInterest)
        
        PoI1=indx(find(cumsum1<PercentsOfInterest(nn)));
        %       figure(2), hold on
        %       plot(F0_SPL_space(2,PoI1),F0_SPL_space(3,PoI1),'k.')
        %       hold off
        PlotPoints=[F0_SPL_space(2,PoI1); F0_SPL_space(3,PoI1)];
        PlotPoints2x=[]; PlotPoints2y=[];
        PlotPointsSize=size(PlotPoints);
        for n1=1:PlotPointsSize(2)
            PlotPoints2x=[PlotPoints2x ...
                PlotPoints(1,n1)-diffF0 PlotPoints(1,n1)+diffF0 ...
                PlotPoints(1,n1)+diffF0 PlotPoints(1,n1)-diffF0];
            PlotPoints2y=[PlotPoints2y ...
                PlotPoints(2,n1)-diffSPL PlotPoints(2,n1)-diffSPL ...
                PlotPoints(2,n1)+diffSPL PlotPoints(2,n1)+diffSPL];
        end % n1
        PlotPoints2x(PlotPoints2x<F0_space(1))=F0_space(1);
        PlotPoints2y(PlotPoints2y<SPL_space(1))=SPL_space(1);
        %     PlotPoints2x=PlotPoints2x(PlotPoints2x>F0_space(1));
        %     PlotPoints2y=PlotPoints2y(PlotPoints2x>F0_space(1));
        %     PlotPoints2y=PlotPoints2y(PlotPoints2y>SPL_space(1));
        %     PlotPoints2x=PlotPoints2x(PlotPoints2y>SPL_space(1));
        
%         if ~isempty(PlotPoints2x)
%             figure(fignum+3)
%             k2 = convhull(PlotPoints2x,PlotPoints2y);
%             if PlotPointsSize(2)>1
%                 try
%                     k = convhull(PlotPoints(1,:),PlotPoints(2,:));
%                     plot(PlotPoints(1,k),PlotPoints(2,k),'b-',...
%                         PlotPoints2x(k2),PlotPoints2y(k2),'r:')
%                     tmp=not(PlotPoints(1,k)==F0_space(1)&PlotPoints(2,k)==SPL_space(1));
%                     b=fit_ellipse(PlotPoints(1,k(tmp)),PlotPoints(2,k(tmp)));
%                     try
%                         if isempty(b.a)
%                             plot(PlotPoints(1,:),PlotPoints(2,:),'b-',...
%                                 PlotPoints2x(k2),PlotPoints2y(k2),'r:')
%                             k=k2;
%                             tmp=not(PlotPoints2x(k2)==F0_space(1)&PlotPoints2y(k2)==SPL_space(1));
%                             b=fit_ellipse(PlotPoints2x(k2(tmp)),PlotPoints2y(k2(tmp)));
%                         end
%                     catch
%                         if isempty(b)
%                             plot(PlotPoints(1,:),PlotPoints(2,:),'b-',...
%                                 PlotPoints2x(k2),PlotPoints2y(k2),'r:')
%                             k=k2;
%                             tmp=not(PlotPoints2x(k2)==F0_space(1)&PlotPoints2y(k2)==SPL_space(1));
%                             b=fit_ellipse(PlotPoints2x(k2(tmp)),PlotPoints2y(k2(tmp)));
%                         end
%                     end
%                 catch
%                     plot(PlotPoints(1,:),PlotPoints(2,:),'b-',...
%                         PlotPoints2x(k2),PlotPoints2y(k2),'r:')
%                     k=k2;
%                     tmp=not(PlotPoints2x(k2)==F0_space(1)&PlotPoints2y(k2)==SPL_space(1));
%                     b=fit_ellipse(PlotPoints2x(k2(tmp)),PlotPoints2y(k2(tmp)));
%                 end
%                 
%             else
%                 plot(PlotPoints(1,:),PlotPoints(2,:),'b-',...
%                     PlotPoints2x(k2),PlotPoints2y(k2),'r:')
%                 k=k2;
%                 tmp=not(PlotPoints2x(k2)==F0_space(1)&PlotPoints2y(k2)==SPL_space(1));
%                 b=fit_ellipse(PlotPoints2x(k2(tmp)),PlotPoints2y(k2(tmp)));
%             end
%             
%             tmp=PercentsOfInterest(nn);
%             tmp2=round(10*length(PlotPoints)*(4*diffSPL*diffF0))/10;
%             text(0.995*max(PlotPoints(1,:)),0.99*max(PlotPoints(2,:)),...
%                 [' \leftarrow ' num2str(tmp) '%, A=' num2str(tmp2) 'stSPL'])%,'FontSize',18)
%             
            
%             % the ellipse
%             if PlotPointsSize(2)>1&&isempty(b)%~isempty(b.a)
%                 R = [ cos(b.phi) sin(b.phi); -sin(b.phi) cos(b.phi) ];
%                 theta_r         = linspace(0,2*pi);
%                 ellipse_x_r     = b.X0 + b.a*cos( theta_r );
%                 ellipse_y_r     = b.Y0 + b.b*sin( theta_r );
%                 rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];
%                 a=fit_ellipse(rotated_ellipse(1,:),rotated_ellipse(2,:));
%                 
%                 figure(fignum+4),
%                 plot(rotated_ellipse(1,:),rotated_ellipse(2,:),'r',...
%                     b.X0_in,b.Y0_in,'ko')
%                 
%                 xlabel(xlabtext); ylabel(ylabtext);
%                 
%                 
%                 % trim elips if beyond the Fo and SPL space then use Polyarea to get area
%                 if sum(rotated_ellipse(1,:)<F0_space(1))>0&&sum(rotated_ellipse(2,:)<SPL_space(1))>0
%                     tmp1=rotated_ellipse(1,:);
%                     tmp1(tmp1<F0_space(1))=F0_space(1);
%                     tmp2=rotated_ellipse(2,:);
%                     tmp2(tmp2<SPL_space(1))=SPL_space(1);
%                     tmp=polyarea(tmp1,tmp2);
%                     tmp=round(tmp); %area
%                 else
%                     tmp=round(pi*b.a*b.b); %area
%                 end
%                 tmp3=round(100*(1/tan(-b.phi)))/100;
%                 
%                 tmp2=round(length(PlotPoints)*(4*diffSPL*diffF0));
%                 
%                 text(max(rotated_ellipse(1,:)),max(rotated_ellipse(2,:)),...
%                     ['A_{el}=' num2str(tmp) ...
%                     ', A=' num2str(tmp2) ' stSPL' ...
%                     '; slp=' num2str(tmp3) ' dB/st'],'FontSize',8)              
%             else
%                 tmp=0; %area
%                 tmp3=0;
%             end
%         end
    end % nn
    %%
%     figure(fignum+3)
%     plot(F0_SPL_space(2,indx(1)),F0_SPL_space(3,indx(1)),'g.')
%     hold off
%     axis([F0_space(1),F0_space(end),SPL_space(1),SPL_space(end)])
%     xlabel(xlabtext); ylabel(ylabtext);
%     tmp=num2str(round(10*F0_SPL_space(2,indx(1)))/10);
%     tmp2=num2str(round(10*F0_SPL_space(3,indx(1))/10));
%     title( [strrep(filename, '_', '-') ': mode= ' tmp ' st, ' tmp2 ' dB SPL'])
%     
%     figure(fignum+4)
%     plot(F0_SPL_space(2,indx(1)),F0_SPL_space(3,indx(1)),'g.')
%     hold off
%     axis([F0_space(1),F0_space(end),SPL_space(1),SPL_space(end)])
%     xlabel(xlabtext); ylabel(ylabtext);
%     tmp=num2str(round(10*F0_SPL_space(2,indx(1)))/10);
%     tmp2=num2str(round(10*F0_SPL_space(3,indx(1))/10));
%     title( [strrep(filename, '_', '-') ': mode= ' tmp ' st, ' tmp2 ' dB SPL'])
    
    




    
    %   figure
    %   plot(rotated_ellipse(1,:),rotated_ellipse(2,:),':',tmp1,tmp2)
    %   polyarea(rotated_ellipse(1,:),rotated_ellipse(2,:))
    %   polyarea(tmp1,tmp2)
    
    %   %% fit ellipse
    %     b=fit_ellipse(PlotPoints(1,k),PlotPoints(2,k));
    %     % the ellipse
    %     R = [ cos(b.phi) sin(b.phi); -sin(b.phi) cos(b.phi) ];
    %     theta_r         = linspace(0,2*pi);
    %     ellipse_x_r     = b.X0 + b.a*cos( theta_r );
    %     ellipse_y_r     = b.Y0 + b.b*sin( theta_r );
    %     rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];
    %     a=fitellipse(rotated_ellipse(1,:),rotated_ellipse(2,:));
    %
    % %     hypot           = (linspace(0,b.long_axis)-b.long_axis/2);
    % %     major_x         = hypot*cos(a(5))+b.X0_in;
    % %     major_y         = hypot*sin(a(5))+b.Y0_in;
    % %     hypot           = (linspace(0,b.short_axis)-b.short_axis/2);
    % %     minor_x         = hypot*cos(a(5))+b.X0_in;
    % %     minor_y         = hypot*sin(a(5))+b.Y0_in;
    %
    %     figure(4),
    %     plot(...
    %       rotated_ellipse(1,:),rotated_ellipse(2,:),'r',...
    %       F0_SPL_space(2,indx(1)),F0_SPL_space(3,indx(1)),'ks',...
    %       b.X0_in,b.Y0_in,'ko',...
    %       PlotPoints(1,:),PlotPoints(2,:),'b.',...
    %       PlotPoints2x,PlotPoints2y,'g.')
    % %       major_x,major_y,'m',minor_x,minor_y,'m')
    %
    %     %PlotPoints(1,k),PlotPoints(2,k),'b-',...
    %     %PlotPoints2x(k2),PlotPoints2y(k2),'r:',...
    %
    %     xlabel(xlabtext); ylabel(ylabtext);
    %     tmp=round(10*pi*b.a*b.b)/10; %area
    %     tmp2=round(10*length(PlotPoints)*(4*diffSPL*diffF0))/10;
    %     tmp3=round(1000*tan(-b.phi))/1000;
    %     title(['A_{ellipse}=' num2str(tmp) ...
    %       ', A_{calc}=' num2str(tmp2) ' stSPL' ...
    %       '; slope=' num2str(tmp3) ' dB/st'])
    
    %   axis([F0_space(1),600,SPL_space(1),80])
    
    %% save figures and data
    filename=strrep(filename, ' ', '-');
    
%     figure(fignum+4)
%     saveas(gcf,[filename '_VUP_ellipse_areas.jpg'],'jpg')
%     saveas(gcf,[filename '_VUP_ellipse_areas.fig'],'fig')
%     
%     figure(fignum+3)
%     saveas(gcf,[filename '_VUP_areas.jpg'],'jpg')
%     saveas(gcf,[filename '_VUP_areas.fig'],'fig')
    
    figure(fignum+1)
    saveas(gcf,[filename '_VUP_3D.fig'],'fig')
    saveas(gcf,[filename '_VUP_3D.jpg'],'jpg')
    
    figure(fignum+2)
    %   axis([F0_space(1),600,SPL_space(1),80])
    saveas(gcf,[filename '_VUP_contour.jpg'],'jpg')
    saveas(gcf,[filename '_VUP_contour.fig'],'fig')
    
    
    save([filename '_VUP.mat'])
    
    disp( '********  VUP complete  *************')
%     disp(['  VUP complete  '])
%     disp( '*********************')
    
else
%     disp('.')
    disp( '**********  error with file!!  ***********')
%     disp(['  error with file!!  '])
%     disp( '*********************')
%     disp('.')
    figure(fignum+1), clf
    title( [strrep(filename, '_', '-') ': ERROR WITH FILE!!!!!'])
    saveas(gcf,[filename(1:end-4) '_VUP_areas.jpg'],'jpg')
    saveas(gcf,[filename(1:end-4) '_VUP_areas.fig'],'fig')
end



%% find where each individual contour is.
%   contour_strt=[];
%   base1=C(1,1);
%   n=0;
%   tmp=find(round(100*C(1,1:end))==round(100*(n+1)*base1));
%   contour_strt(1:length(tmp),:)=[tmp; ones(1,length(tmp)).*(n+1)*base1]';
%   n=n+1;
%   quit1=1;
%   while quit1
%     tmp=find(round(100*C(1,contour_strt(end,1):end))==round(100*(n+1)*base1));
%     if ~isempty(tmp>0)
%       contour_strt=...
%         [contour_strt; [tmp+contour_strt(end,1)-1; ones(1,length(tmp)).*(n+1)*base1]'];
%       n=n+1;
%     else
%       quit1=0;
%     end
%   end
%   %%% add the final marker
%   contour_strt=...
%     [contour_strt; [length(C)+1 contour_strt(end,2)]];
%
%   %% find number of occurances in a contoured area
%   clear A
%   Tot_Occur=sum(F0_SPL_mat(:,3)); % total occurances
%   Max_Occur=max(F0_SPL_mat(:,3)); % Max occurances
%   figure(2), clf, hold on
%   %%%%%%%%%%%%%%%%% find occurances within each individual contour
%   for n=1:length(contour_strt)-1
%     tmpi=contour_strt(n:n+1,1);
%     fo_x=C(1,tmpi(1)+1:tmpi(2)-1);
%     spl_y=C(2,tmpi(1)+1:tmpi(2)-1);
%     if fo_x(1)~=fo_x(end) || spl_y(1)~=spl_y(end)
%       fo_x=[fo_x fo_x(1)];
%       spl_y=[spl_y spl_y(1)];
%     end
%     %       plot(fo_x,spl_y)
%     [in on] = inpolygon(F0_SPL_mat(:,1),F0_SPL_mat(:,2),fo_x,spl_y);
%     C_TotalOccurance(n) = sum(F0_SPL_mat(in,3));
%     C_PerOfTotOccur(n)  = C_TotalOccurance(n)/Tot_Occur;
%     C_Area(n)           = polyarea(fo_x,spl_y);
%     C_Level(n)          = contour_strt(n,2);
%     C_PerOfMaxLevel(n)  = contour_strt(n,2)/F0_SPL_mat(F0_SPL_mat_maxi(3),3);
%     C_Averages(n,:)     = [mean(F0_SPL_mat(in,1)) mean(F0_SPL_mat(in,2))];
%     C_AvgReal(n,:)      = [sum(F0_SPL_mat(in,1).*F0_SPL_mat(in,3))/sum(F0_SPL_mat(in,3)) ...
%       sum(F0_SPL_mat(in,2).*F0_SPL_mat(in,3))/sum(F0_SPL_mat(in,3))];
%   end
%   hold off
%   %% find contour areas of interest
%   PercentsOfInterest=[0.682, 0.5 0.25 0.1];
%   ContoursOfInterest=[];
%   err0=0.00001;err=err0;
%   n=0;nogo=0;
%   for m=1:length(PercentsOfInterest)
%     while nogo==0
%       tmp=find(C_PerOfTotOccur>(PercentsOfInterest(m)-err)&...
%         C_PerOfTotOccur<(PercentsOfInterest(m)+err));
%       if length(tmp)<1
%         err=err+err0;
%         nogo=0;
%       else
%         nogo=1;
%         ContoursOfInterest=[ContoursOfInterest...
%           [tmp(1) ; 0.*tmp(1)+PercentsOfInterest(m)]];
%         err=err0;
%       end
%     end
%     nogo=0;
%   end
%
%   %% plot each contour of interest
%   figure(2), clf, hold on
%   %%%%%%%%%%%%%%%%% find where each individual contour is.
%   for n=1:length(ContoursOfInterest(1,:))
%     nn=(ContoursOfInterest(1,n));
%     tmpi=contour_strt(nn:nn+1,1);
%     fo_x=C(1,tmpi(1)+1:tmpi(2)-1);
%     spl_y=C(2,tmpi(1)+1:tmpi(2)-1);
%     if fo_x(1)~=fo_x(end) || spl_y(1)~=spl_y(end)
%       fo_x=[fo_x fo_x(1)];
%       spl_y=[spl_y spl_y(1)];
%     end
%     plot(fo_x,spl_y)
%     [fo_max fo_maxi]=max(fo_x);
%     tmp=round(10*100*C_PerOfTotOccur(nn))/10;
%     tmp2=round(10*100*C_Area(nn))/10;
%     plot(C_Averages(nn,1),C_Averages(nn,2),'g.')
%     plot(C_AvgReal(nn,1),C_AvgReal(nn,2),'k.')
%     text(fo_x(fo_maxi),spl_y(fo_maxi),...
%       [' \leftarrow ' num2str(tmp) '%, A=' num2str(tmp2) 'stSPL'])%,'FontSize',18)
%     pause(.11)
%     %     keyboard
%   end
%   plot(F0_SPL_mat(F0_SPL_mat_maxi(3),1),F0_SPL_mat(F0_SPL_mat_maxi(3),2),'r+')
%   plot(sum(F0_SPL_mat(:,1).*F0_SPL_mat(:,3))/sum(F0_SPL_mat(:,3)),...
%     sum(F0_SPL_mat(:,2).*F0_SPL_mat(:,3))/sum(F0_SPL_mat(:,3)),'k+')
%   plot(fo_x,spl_y,'r:')
%   hold off
%   xlabel(xlabtext); ylabel(ylabtext);
%   title( [strrep(filename, '_', '-') ': 10% Area= ' num2str(tmp2) ' stSPL'])
%   %% save figures and data
%   saveas(gcf,[filename(1:end-4) '_VUP_areas.jpg'],'jpg')
%   saveas(gcf,[filename(1:end-4) '_VUP_areas.fig'],'fig')
%   save([filename '_VUP.mat'])
%
% else
%   disp('.')
%   disp( '*********************')
%   disp(['  error with file!!  '])
%   disp( '*********************')
%   disp('.')
%   figure(111), clf
%   title( [strrep(filename, '_', '-') ': ERROR WITH FILE!!!!!'])
%   saveas(gcf,[filename(1:end-4) '_VUP_areas.jpg'],'jpg')
%   saveas(gcf,[filename(1:end-4) '_VUP_areas.fig'],'fig')
% end





%%
% % *************************************************************************
% % Creates an array for use in surface plots
% % *************************************************************************
% n = 1;
% for m = 1:numel(SPL_array)
%     if SPL_array(m) == F0_SPL_mat(n,2)
%         F0_SPL_array(m) = F0_SPL_mat(n,3);
%         n = n + 1;
%     end
%     if n > length(F0_SPL_mat) break; end
% end
%
% % *************************************************************************
% % Set plotting parameters
% % *************************************************************************
%
% if n>1
%
% F0_axes = [min(F0_space) max(F0_space)];
% SPL_axes = [min(SPL_voiced_sorted) SPL_space(end)];
% view_angle = [40 35];
%
% plot(F0_SPL_array)
%
%
% F0_SPL_array0=F0_SPL_array;
% N=6;
% slp=mean(diff(F0_SPL_array(N:N+6,:)));
% for n=1:N
% % slp=(F0_SPL_array(N,:)-F0_SPL_array(1,:))/(N-1);
% F0_SPL_array0(n,:)=slp.*(n-1)+F0_SPL_array(N-1,:)-(N)*slp;
% end
% plot(F0_SPL_array0)
%
%
% figure(2)
% h = surf(F0_space,SPL_space,F0_SPL_array0,'EdgeColor','none',...
%    'FaceLighting','phong','FaceColor','interp','DiffuseStrength',1,...
%     'SpecularStrength',0.3);
% xlim(F0_axes+[0 -600])
% ylim(SPL_axes+[0 -10])
% xlabel(xlabtext); ylabel(ylabtext); zlabel('# Occurrences');
% grid on
% view(view_angle)
%
% x = linspace(0,1,64);
% cmap = colormap(winter);
% colormap(cmap(64:-1:1,:));
% alpha('z');
% alpha_map = x.^0.1; % this removes most of the low count occurrences
% alphamap(alpha_map);
% lightangle(146,20)
% title(['# occur = ' num2str(round(sum(sum(F0_SPL_array0))))])
% saveas(gcf,[filename(1:end-4) '_3D'],'fig')
%
% save VUP


%%  Plotting stuff


% figure(111)
% % *************************************************************************
% % Plots the data as x = F0, y = SPL, and z = occurrences
% % *************************************************************************
% subplot 221
% h = scatter3(F0_SPL_mat(:,1),F0_SPL_mat(:,2),F0_SPL_mat(:,3),10,[0 .9 .2]);
% xlabel(xlabtext); ylabel(ylabtext); zlabel('# Occurrences');
% xlim(F0_axes)
% ylim(SPL_axes)
% % get(gca,'XTick')
% grid on
% view(view_angle)
%
% % *************************************************************************
% % Plots the matrix F0_SPL_array as a surface
% % *************************************************************************
% % figure
% subplot 222
% h = surf(F0_space,SPL_space,F0_SPL_array,'EdgeColor','none',...
%    'FaceLighting','phong','FaceColor','interp','DiffuseStrength',1,...
%     'SpecularStrength',0.3);
% xlim(F0_axes+[0 -600])
% ylim(SPL_axes+[0 -10])
% xlabel(xlabtext); ylabel(ylabtext); zlabel('# Occurrences');
% grid on
% view(view_angle)
%
%
% % Sets the color for the data and adjusts lighting and transparency
% x = linspace(0,1,64);
% % color_map = x.^0.01;
% cmap = colormap(winter);
% % cmap(:,2) = color_map';
% % colormap(cmap)
% colormap(cmap(64:-1:1,:));
% alpha('z');
% % x = 0:1/64:1;
% alpha_map = x.^0.1; % this removes most of the low count occurrences
% % plot(alpha_map) % shows how the transparency is applied
% alphamap(alpha_map);
% lightangle(146,20)
% % light('Position',[10 0.3 70],'Style','infinite');
%
%
% % *************************************************************************
% % Plots the matrix F0_SPL_array as a 3-D contour
% % *************************************************************************
% subplot 223
% contour3(F0_space,SPL_space,F0_SPL_array,60)
% xlim(F0_axes+[0 -600])
% ylim(SPL_axes+[0 -10])
% xlabel(xlabtext); ylabel(ylabtext); zlabel('# Occurrences');
% grid on
% view(view_angle)
%
%
% % *************************************************************************
% % Plots the matrix F0_SPL_array as a 2-D filled contour
% % *************************************************************************
% subplot 224
% contourf(F0_space,SPL_space,F0_SPL_array,20)
% xlim(F0_axes+[0 -600])
% ylim(SPL_axes+[0 -10])
% xlabel(xlabtext); ylabel(ylabtext); zlabel('# Occurrences');
% colorbar
% % title(['                           '...
% %        '                           '...
% %        '                           '...
% %        'Occurrences'])
%
% %%



function ellipse_t = fit_ellipse( x,y,axis_handle )
%
% fit_ellipse - finds the best fit to an ellipse for the given set of points.
%
% Format:   ellipse_t = fit_ellipse( x,y,axis_handle )
%
% Input:    x,y         - a set of points in 2 column vectors. AT LEAST 5 points are needed !
%           axis_handle - optional. a handle to an axis, at which the estimated ellipse
%                         will be drawn along with it's axes
%
% Output:   ellipse_t - structure that defines the best fit to an ellipse
%                       a           - sub axis (radius) of the X axis of the non-tilt ellipse
%                       b           - sub axis (radius) of the Y axis of the non-tilt ellipse
%                       phi         - orientation in radians of the ellipse (tilt)
%                       X0          - center at the X axis of the non-tilt ellipse
%                       Y0          - center at the Y axis of the non-tilt ellipse
%                       X0_in       - center at the X axis of the tilted ellipse
%                       Y0_in       - center at the Y axis of the tilted ellipse
%                       long_axis   - size of the long axis of the ellipse
%                       short_axis  - size of the short axis of the ellipse
%                       status      - status of detection of an ellipse
%
% Note:     if an ellipse was not detected (but a parabola or hyperbola), then
%           an empty structure is returned

% =====================================================================================
%                  Ellipse Fit using Least Squares criterion
% =====================================================================================
% We will try to fit the best ellipse to the given measurements. the mathematical
% representation of use will be the CONIC Equation of the Ellipse which is:
%
%    Ellipse = a*x^2 + b*x*y + c*y^2 + d*x + e*y + f = 0
%
% The fit-estimation method of use is the Least Squares method (without any weights)
% The estimator is extracted from the following equations:
%
%    g(x,y;A) := a*x^2 + b*x*y + c*y^2 + d*x + e*y = f
%
%    where:
%       A   - is the vector of parameters to be estimated (a,b,c,d,e)
%       x,y - is a single measurement
%
% We will define the cost function to be:
%
%   Cost(A) := (g_c(x_c,y_c;A)-f_c)'*(g_c(x_c,y_c;A)-f_c)
%            = (X*A+f_c)'*(X*A+f_c)
%            = A'*X'*X*A + 2*f_c'*X*A + N*f^2
%
%   where:
%       g_c(x_c,y_c;A) - vector function of ALL the measurements
%                        each element of g_c() is g(x,y;A)
%       X              - a matrix of the form: [x_c.^2, x_c.*y_c, y_c.^2, x_c, y_c ]
%       f_c            - is actually defined as ones(length(f),1)*f
%
% Derivation of the Cost function with respect to the vector of parameters "A" yields:
%
%   A'*X'*X = -f_c'*X = -f*ones(1,length(f_c))*X = -f*sum(X)
%
% Which yields the estimator:
%
%       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%       |  A_least_squares = -f*sum(X)/(X'*X) ->(normalize by -f) = sum(X)/(X'*X)  |
%       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% (We will normalize the variables by (-f) since "f" is unknown and can be accounted for later on)
%
% NOW, all that is left to do is to extract the parameters from the Conic Equation.
% We will deal the vector A into the variables: (A,B,C,D,E) and assume F = -1;
%
%    Recall the conic representation of an ellipse:
%
%       A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0
%
% We will check if the ellipse has a tilt (=orientation). The orientation is present
% if the coefficient of the term "x*y" is not zero. If so, we first need to remove the
% tilt of the ellipse.
%
% If the parameter "B" is not equal to zero, then we have an orientation (tilt) to the ellipse.
% we will remove the tilt of the ellipse so as to remain with a conic representation of an
% ellipse without a tilt, for which the math is more simple:
%
% Non tilt conic rep.:  A`*x^2 + C`*y^2 + D`*x + E`*y + F` = 0
%
% We will remove the orientation using the following substitution:
%
%   Replace x with cx+sy and y with -sx+cy such that the conic representation is:
%
%   A(cx+sy)^2 + B(cx+sy)(-sx+cy) + C(-sx+cy)^2 + D(cx+sy) + E(-sx+cy) + F = 0
%
%   where:      c = cos(phi)    ,   s = sin(phi)
%
%   and simplify...
%
%       x^2(A*c^2 - Bcs + Cs^2) + xy(2A*cs +(c^2-s^2)B -2Ccs) + ...
%           y^2(As^2 + Bcs + Cc^2) + x(Dc-Es) + y(Ds+Ec) + F = 0
%
%   The orientation is easily found by the condition of (B_new=0) which results in:
%
%   2A*cs +(c^2-s^2)B -2Ccs = 0  ==> phi = 1/2 * atan( b/(c-a) )
%
%   Now the constants   c=cos(phi)  and  s=sin(phi)  can be found, and from them
%   all the other constants A`,C`,D`,E` can be found.
%
%   A` = A*c^2 - B*c*s + C*s^2                  D` = D*c-E*s
%   B` = 2*A*c*s +(c^2-s^2)*B -2*C*c*s = 0      E` = D*s+E*c
%   C` = A*s^2 + B*c*s + C*c^2
%
% Next, we want the representation of the non-tilted ellipse to be as:
%
%       Ellipse = ( (X-X0)/a )^2 + ( (Y-Y0)/b )^2 = 1
%
%       where:  (X0,Y0) is the center of the ellipse
%               a,b     are the ellipse "radiuses" (or sub-axis)
%
% Using a square completion method we will define:
%
%       F`` = -F` + (D`^2)/(4*A`) + (E`^2)/(4*C`)
%
%       Such that:    a`*(X-X0)^2 = A`(X^2 + X*D`/A` + (D`/(2*A`))^2 )
%                     c`*(Y-Y0)^2 = C`(Y^2 + Y*E`/C` + (E`/(2*C`))^2 )
%
%       which yields the transformations:
%
%           X0  =   -D`/(2*A`)
%           Y0  =   -E`/(2*C`)
%           a   =   sqrt( abs( F``/A` ) )
%           b   =   sqrt( abs( F``/C` ) )
%
% And finally we can define the remaining parameters:
%
%   long_axis   = 2 * max( a,b )
%   short_axis  = 2 * min( a,b )
%   Orientation = phi
%
%

% initialize
orientation_tolerance = 1e-3;

% empty warning stack
warning( '' );

% prepare vectors, must be column vectors
x = x(:);
y = y(:);

% remove bias of the ellipse - to make matrix inversion more accurate. (will be added later on).
mean_x = mean(x);
mean_y = mean(y);
x = x-mean_x;
y = y-mean_y;

% the estimation for the conic equation of the ellipse
X = [x.^2, x.*y, y.^2, x, y ];
a = sum(X)/(X'*X);

% check for warnings
if ~isempty( lastwarn )
    disp( 'stopped because of a warning regarding matrix inversion' );
    ellipse_t = [];
    return
end

% extract parameters from the conic equation
[a,b,c,d,e] = deal( a(1),a(2),a(3),a(4),a(5));

% remove the orientation from the ellipse
if ( min(abs(b/a),abs(b/c)) > orientation_tolerance )
    
    orientation_rad = 1/2 * atan( b/(c-a) );
    cos_phi = cos( orientation_rad );
    sin_phi = sin( orientation_rad );
    [a,b,c,d,e] = deal(...
        a*cos_phi^2 - b*cos_phi*sin_phi + c*sin_phi^2,...
        0,...
        a*sin_phi^2 + b*cos_phi*sin_phi + c*cos_phi^2,...
        d*cos_phi - e*sin_phi,...
        d*sin_phi + e*cos_phi );
    [mean_x,mean_y] = deal( ...
        cos_phi*mean_x - sin_phi*mean_y,...
        sin_phi*mean_x + cos_phi*mean_y );
else
    orientation_rad = 0;
    cos_phi = cos( orientation_rad );
    sin_phi = sin( orientation_rad );
end

% check if conic equation represents an ellipse
test = a*c;
switch (1)
    case (test>0),  status = '';
    case (test==0), status = 'Parabola found';  warning( 'fit_ellipse: Did not locate an ellipse' );
    case (test<0),  status = 'Hyperbola found'; warning( 'fit_ellipse: Did not locate an ellipse' );
end

% if we found an ellipse return it's data
if (test>0)
    
    % make sure coefficients are positive as required
    if (a<0), [a,c,d,e] = deal( -a,-c,-d,-e ); end
    
    % final ellipse parameters
    X0          = mean_x - d/2/a;
    Y0          = mean_y - e/2/c;
    F           = 1 + (d^2)/(4*a) + (e^2)/(4*c);
    [a,b]       = deal( sqrt( F/a ),sqrt( F/c ) );
    long_axis   = 2*max(a,b);
    short_axis  = 2*min(a,b);
    
    % rotate the axes backwards to find the center point of the original TILTED ellipse
    R           = [ cos_phi sin_phi; -sin_phi cos_phi ];
    P_in        = R * [X0;Y0];
    X0_in       = P_in(1);
    Y0_in       = P_in(2);
    
    % pack ellipse into a structure
    ellipse_t = struct( ...
        'a',a,...
        'b',b,...
        'phi',orientation_rad,...
        'X0',X0,...
        'Y0',Y0,...
        'X0_in',X0_in,...
        'Y0_in',Y0_in,...
        'long_axis',long_axis,...
        'short_axis',short_axis,...
        'status','' );
else
    % report an empty structure
    ellipse_t = struct( ...
        'a',[],...
        'b',[],...
        'phi',[],...
        'X0',[],...
        'Y0',[],...
        'X0_in',[],...
        'Y0_in',[],...
        'long_axis',[],...
        'short_axis',[],...
        'status',status );
end

% check if we need to plot an ellipse with it's axes.
if (nargin>2) & ~isempty( axis_handle ) & (test>0)
    
    % rotation matrix to rotate the axes with respect to an angle phi
    R = [ cos_phi sin_phi; -sin_phi cos_phi ];
    
    % the axes
    %     ver_line        = [ [X0 X0]; Y0+b*[-1 1] ];
    %     horz_line       = [ X0+a*[-1 1]; [Y0 Y0] ];
    %     new_ver_line    = R*ver_line;
    %     new_horz_line   = R*horz_line;
    
    % the ellipse
    theta_r         = linspace(0,2*pi);
    ellipse_x_r     = X0 + a*cos( theta_r );
    ellipse_y_r     = Y0 + b*sin( theta_r );
    rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];
    
    % draw
    hold_state = get( axis_handle,'NextPlot' );
    set( axis_handle,'NextPlot','add' );
    %     plot( new_ver_line(1,:),new_ver_line(2,:),'r' );
    %     plot( new_horz_line(1,:),new_horz_line(2,:),'r' );
    plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );
    set( axis_handle,'NextPlot',hold_state );
end


function CaptureFigVid(ViewZ, FileName,OptionZ)
% CaptureFigVid(ViewZ, FileName,OptionZ)
% Captures a video of the 3D plot in the current axis as it rotates based
% on ViewZ and saves it as 'FileName.mpg'. Option can be specified.
%
% ViewZ:     N-rows with 2 columns, each row are the view angles in
%            degrees, First column is azimuth (pan), Second is elevation
%            (tilt) values outside of 0-360 wrap without error,
%            *If a duration is specified, angles are used as nodes and
%            views are equally spaced between them (other interpolation
%            could be implemented, if someone feels so ambitious).
%            *If only an initial and final view is given, and no duration,
%            then the default is 100 frames.
% FileName:  Name of the file of the produced animation. Because I wrote
%            the program, I get to pick my default of mpg-4, and the file
%            extension .mpg will be appended, even if the filename includes
%            another file extension. File is saved in the working
%            directory.
% (OptionZ): Optional input to specify parameters. The ones I use are given
%            below. Feel free to add your own. Any or all fields can be
%            used
% OptionZ.FrameRate: Specify the frame rate of the final video (e.g. 30;)
% OptionZ.Duration: Specify the length of video in seconds (overrides
%    spacing of view angles) (e.g. 3.5;)
% OptionZ.Periodic: Logical to indicate if the video should be periodic.
%    Using this removed the final view so that when the video repeats the
%    initial and final view are not the same. Avoids having to find the
%    interval between view angles. (e.g. true;)
%
% % % % Example (shown in published results, video attached) % % % %
% figure(171);clf;
% surf(peaks,'EdgeColor','none','FaceColor','interp','FaceLighting','phong')
% daspect([1,1,.3]);axis tight;
% OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
% CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10],'WellMadeVid',OptionZ)
%
% Known issues: MPEG-4 video option only available on Windows machines. See
% fix where the VideoWriter is called.
%
% Getframe is used to capture image and current figure must be on monitor 1
% if multiple displays are used. Does not work if no display is used.
%
% Active windows that overlay the figure are captured in the movie.  Set up
% the current figure prior to calling the function. If you don't specify
% properties, such as tick marks and aspect ratios, they will likely change
% with the rotation for an undesirable effect.

% Cheers, Dr. Alan Jennings, Research assistant professor,
% Department of Aeronautics and Astronautics, Air Force Institute of Technology

%% preliminaries

% initialize optional argument
if nargin<3;     OptionZ=struct([]); end

% check orientation of ViewZ, should be two columns and >=2 rows
if size(ViewZ,2)>size(ViewZ,1); ViewZ=ViewZ.'; end
if size(ViewZ,2)>2
    warning('AJennings:VidWrite',...
        'Views should have n rows and only 2 columns. Deleting extraneous input.');
    ViewZ=ViewZ(:,1:2); %remove any extra columns
end

% Create video object
daObj=VideoWriter(FileName,'MPEG-4'); %my preferred format
% daObj=VideoWriter(FileName); %for default video format.
% MPEG-4 CANNOT BE USED ON UNIX MACHINES
% set values:
% Frame rate
if isfield(OptionZ,'FrameRate')
    daObj.FrameRate=OptionZ.FrameRate;
end
% Durration (if frame rate not set, based on default)
if isfield(OptionZ,'Duration') %space out view angles
    temp_n=round(OptionZ.Duration*daObj.FrameRate); % number frames
    temp_p=(temp_n-1)/(size(ViewZ,1)-1); % length of each interval
    ViewZ_new=zeros(temp_n,2);
    % space view angles, if needed
    for inis=1:(size(ViewZ,1)-1)
        ViewZ_new(round(temp_p*(inis-1)+1):round(temp_p*inis+1),:)=...
            [linspace(ViewZ(inis,1),ViewZ(inis+1,1),...
            round(temp_p*inis)-round(temp_p*(inis-1))+1).',...
            linspace(ViewZ(inis,2),ViewZ(inis+1,2),...
            round(temp_p*inis)-round(temp_p*(inis-1))+1).'];
    end
    ViewZ=ViewZ_new;
end
% space view angles, if needed
if length(ViewZ)==2 % only initial and final given
    ViewZ=[linspace(ViewZ(1,1),ViewZ(end,1)).',...
        linspace(ViewZ(1,2),ViewZ(end,2)).'];
end
% Periodicity
if isfield(OptionZ,'Periodic')&&OptionZ.Periodic==true
    ViewZ=ViewZ(1:(end-1),:); %remove last sample
end
% open object, preparatory to making the video
open(daObj);

%% rotate the axis and capture the video
for kathy=1:size(ViewZ,1)
    view(ViewZ(kathy,:)); drawnow;
    writeVideo(daObj,getframe(gcf)); %use figure, since axis changes size based on view
end

%% clean up
close(daObj);



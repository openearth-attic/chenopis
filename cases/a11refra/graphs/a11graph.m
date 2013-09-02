function [] = a11graph(k,ini,log)

%this function reads data from wave model output file (*.tab/*.sp1),
%observation files and plots this data according to layout 1 (4 panels)

%open standard figure
layout01

%number of panels
npanels   = 4;

%filename
dirname   = 'a11refra';
fname     = 'a11ref01';

%frame text from .ini file
inst      = ini.inst;
wmod      = ini.model;

% number of models
nmod = length(wmod);

modstring  = '';
seperation = '; ';
for imod = 1:nmod
  if (imod ~= nmod)
    modstring = [modstring,wmod{imod},seperation];
  else
    modstring = [modstring,wmod{imod}];
  end
end
modcell{1} = modstring;

logo      = ini.logo;
date      = ini.date;
free      = ini.free;
data      = {inst,modcell,free,date,{dirname,fname}};

%read  wave model output
for imod = 1:nmod
   model{imod} = a11rdmod(k,fname,log,imod);
end

%read observations
observ    = a11rdobs(k,fname,log);

%fill frame
%==========================================================================================
xloc      = {[0.05,0.05,0.05] [0.25,0.05,0.05] [0.05] [0.05] [0.25]};   %position text on x-axes
yloc      = {[0.5] [0.7 0.3] [0.8,0.5,0.2]};                            %position text on y-axes

for ii = [ 1 3:5 ]                                       % fill 5 textboxes
    tag   = ['textbox',num2str(ii)];
    btext = data{ii};
    x_loc = xloc{ii};  
    y_loc = yloc{length(btext)};
    f1    = findobj(gcf,'Tag',tag);
    set(gcf,'CurrentAxes',f1);
  
    if (ii == 2)
        text(x_loc(2),y_loc(1),'Wave model(s):');
    end
   
    for jj = 1:length(btext)      
        if(ii == 5)
            text(0.05,0.7,'Case');
            text(0.30,0.7,':');
            text(0.35,0.7,btext(1),'FontWeight','bold');
            text(0.05,0.3,'Fig.');
            text(0.30,0.3,':');
            text(0.35,0.3,btext(2),'FontWeight','bold');
         else
            t = text(x_loc(jj),y_loc(jj),btext(jj)); 
         end
     end
end

%institutes logo
lname = char(logo{1});
a     = exist(lname);

if(a == 2)
    f1    = findobj(gcf,'Tag','logo');
    f2    = imread(lname);
    set(gcf,'CurrentAxes',f1);
    imagesc(f2);
    set(f1,'Visible','off');
end

%draw figures
%==========================================================================================
%line properties: wave model, observations
c          = {[1 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0]};
cb         = {[1 0 0],[0 0 0],[1 0 0]}; %color bottom; bottom, waterlevel, zerolevel
lt         = {'-','-','-','-','-'};
ltb        = {'-','-','--'};            %line type bottom; bottom, waterlevel, zerolevel
width      = [2,0.5,0.5,0.5,0.5];
widthb     = [1.5,1,.3];                %line width bottom; bottom, waterlevel, zerolevel
mark       = {'none','V','^','*','o'};
size       = [6,4,4,4,4];               %default = 6

%x-axes properties
xlim       = [0 4000];
xtick      = [0 500 1000 1500 2000 2500 3000 3500 4000];
xticklabel = [0 500 1000 1500 2000 2500 3000 3500 4000];

%y-axes properties
ylim       = {[-25 5],[0 2.5],[0 12],[80 130]};
ytick      = {[-25:5:5],[0:.5:2.5],[0:2:12],[80:10:130]};
yticklabel = {[-25:5:5],[0:.5:2.5],[0:2:12],[80:10:130]};

%columns to use frome wave model output
modcol    = [1,1,1,1;2,3,4,8]; %first column is the x-axes; second column is the y-axes
obscol    = [1,1,1,1;4,2,1,3]; %first column is the x-axes; second column is the y-axes
nfigs     = 0; %initialise the number of figs
%chart data
%-----------------------------------------------------------------
%wave model
for imod = 1:nmod
  mod = model{imod}.tab; 
  mod0 = mod{1};
  
%number of steps
  nmarks    = 15;    %number of marks
  if (length(mod0(:,:)) > nmarks)
    step = floor(length(mod0(:,:))/nmarks);
  else
    step = 1;
  end

  clear check;
  check                  = find(mod0(:,modcol(2,1)) == -9 | mod0(:,modcol(2,1)) == -99 | ...
                                mod0(:,modcol(2,1)) == -999); 
  mod0(check,modcol(:,:)) = NaN;                   % replace dummy values and set them to zero
  mod0(:,modcol(2,1))     = -mod0(:,modcol(2,1));  % make depth relative to zero
  mod2{1} = mod0;
  mod1{imod} = mod2{1};

end

%observations
obs                    = observ.data;
bot                    = observ.bot;
header                 = model{1}.header;

%set legend
%-----------------------------------------------------------------
%legend text

legtext {1} = 'Analytical solution'; %legend text
for imod = 1:nmod
   legtext {imod+1} = wmod{imod}; 
end
f1      = findobj(gcf,'Tag','legend');                    
set(gcf,'CurrentAxes',f1);
xlegend = [0.  0.  0.  0.5 0.5];             %x position legend text
ylegend = [0.9 0.5 0.1 0.9 0.5];             %y position legend text

for ii = 1:length(legtext)
    xs        = xlegend(ii):0.05:(xlegend(ii)+.1);
    ys        = ylegend(ii).*ones(length(xs));
    l         = strcat('l',num2str(ii));
    l         = line(xs,ys);
    set(l,'Color',c{ii});                            
    set(l,'LineStyle',lt{ii},'LineWidth',width(ii));                    
    set(l,'Marker',mark{ii},'MarkerSize',size(ii));
    text(xlegend(ii)+.15,ylegend(ii),legtext{ii});
end

%plot panels
%-----------------------------------------------------------------
nfigs    =    nfigs + 1;
for ii = 1:npanels
   panel = strcat('panel',num2str(ii));
   f1    = findobj(gcf,'Tag',panel);
   set(gcf,'CurrentAxes',f1);
   set(f1,'Visible','on');                        
   set(f1,'Xlim',xlim);                            
   set(f1,'Xtick',xtick,'XtickLabel',xticklabel);                        
   set(f1,'Ylim',ylim{ii});                    
   set(f1,'Ytick',ytick{ii},'YtickLabel',yticklabel{ii});                        
   for imod = 1:nmod
     mod3 = mod1{imod}; 
     mod  = mod3;
     x     = {mod(:,modcol(1,ii));obs(:,obscol(1,ii))};
     x2    = {mod(1:step:end,modcol(1,ii));obs(:,obscol(1,ii))};
     y     = {mod(:,modcol(2,ii));obs(:,obscol(2,ii))};
     y2    = {mod(1:step:end,modcol(2,ii));obs(:,obscol(2,ii))};
     xx{1}  = x{2};
     xx2{1} = x2{2};
     yy{1}  = y{2};
     yy2{1} = y2{2};
     xx{imod+1}  = x{1};
     xx2{imod+1} = x2{1};
     yy{imod+1}  = y{1};
     yy2{imod+1} = y2{1};
   end   
     
   %ylabel
   %-----------------------------------------------------------------
   labely = header(modcol(2,ii));
   ylabel(labely);
   
   %xlabel
   %-----------------------------------------------------------------
   if ii == 4
      labelx = header(modcol(1,1));
      xlabel(labelx);
   end    
   
   %secundary axes
   %-----------------------------------------------------------------
   ax1 = gca;
   box off;
   ax2 = axes('Position',get(ax1,'Position'),...
     'XaxisLocation','top',...
     'YaxisLocation','right',...
     'Color','none');
   set(ax2,'Xlim',xlim);   
   set(ax2,'Ylim',ylim{ii});
   set(ax2,'Xtick',xlim,'Ytick',ylim{ii});
   set(ax2,'XtickLabel',[],'YtickLabel',[]);
      
   if ii == 1    %plot bottom
      for jj = 1:2
           lh = strcat('lh',num2str(jj));
           lb = strcat('lb',num2str(jj));
           lh = line(bot(:,1),bot(:,jj+1));
           set(lh,'Color',cb{jj});                
           set(lh,'LineStyle',ltb{jj},'LineWidth',widthb(jj)); 
           lb = line(x2{ii},y2{ii});
           set(lb,'Color',c{ii});                
           set(lb,'LineStyle',lt{ii},'LineWidth',width(ii));                    
           set(lb,'MarkerSize',size(ii),'Marker',mark{ii});     
       end
    
   else 
     if (ii == 3) % panel 3; no observations for Tpeak available
       for jj = 2:nmod+1
           lh    = strcat('lh',num2str(jj)); 
           lh = line(xx{jj},yy{jj});
           set(lh,'Color',c{jj});                
           set(lh,'LineStyle',lt{jj},'LineWidth',width(jj));                    
           lh = line(xx2{jj},yy2{jj});
           set(lh,'Color',c{jj});                
           set(lh,'LineStyle','none');                    
           set(lh,'MarkerSize',size(jj),'Marker',mark{jj});
       end   
     else  
       for jj = 1:nmod+1
           lh    = strcat('lh',num2str(jj)); 
           lh = line(xx{jj},yy{jj});
           set(lh,'Color',c{jj});                
           set(lh,'LineStyle',lt{jj},'LineWidth',width(jj));                    
           lh = line(xx2{jj},yy2{jj});
           set(lh,'Color',c{jj});                
           set(lh,'LineStyle','none');                    
           set(lh,'MarkerSize',size(jj),'Marker',mark{jj});
       end   
     end
   end     
end

%save figures
savefig(ini,dirname,nfigs,k)

function [observ] = a11rdobs(k,fname,log)

%cross section plots
%read observation files

%extention observation files
ext     = {'.ana'};

%number of columns in observation files
num_col = [7];

%total number of headerlines in observation files
tot_lin = [5];

%check platform UNIX or DOS
if(k)
   str = ['cd ../observ'];   
else
   str = ['cd ..\observ'];
end

for ii = 1:1
   eval(str);
   f1name    = strcat(fname,ext{ii});
   
   % check if file exist
    a = exist(f1name);
   
   if(a == 2)
       fid = fopen(f1name,'r');
      
       for jj = 1:tot_lin(ii)
           line = fgetl(fid);      
       end
       
       [obs]    = fscanf(fid,'%f',[num_col(ii) inf]);
       obs      = obs';
       fclose(fid);
    else
       fprintf(log,'%s %s %s\n','Observation file ',f1name,' does not exist');
       obs(:,1:num_col(ii)) = NaN;
   end   
end

%read bottom file
%============================================================================
%startpoint bottom
x0      = 0;

%endpoint bottom
xn      = 4000;

%number of meshes
minp    = 8;

%mesh size
dx      = 500;

options = {'linear','nonlinear'};
thisrun = options{1};
switch(thisrun)

case('linear')
   %water level
   mwl = 0;
   
   %depth at x0
   depth0 = -20;
   
   %depth at xn
   depthn = 0;
   
   %zero level
   level0     = 0; 
   [bottom] = [x0 depth0 mwl level0;xn depthn mwl level0];   
   
case('nonlinear')
   
   if(k)      
       str =['cd ../data_in'];   
    else
       str =['cd ..\data_in'];
   end
   eval(str);
   
   %number of headerlines in .bot file
   nlbot   = 1;
   
   %number of headerlines in .bot file
   ncbot   = 1;

   fid = fopen('a11refr.bot');
   for ii = 1:nlbot
      line = fgetl(fid);
   end
   
   bot     = fscanf(fid,'%f',[ncbot inf]);
   bot     = -bot';
   fclose(fid)
   
   %number of headerlines in .lev file
   nllev   = 1;
   
   %number of headerlines in .lev file
   nclev   = 1;
   
   levname = 'a11refra.lev';
   
   % check if file exist
    a = exist(f1name);
   
   %zerolevel
   
   level0 = 0;
   
   %read water level file
   if(a == 2)
       fid = fopen(levname);
       for ii = 1:nllev
           line = fgetl(fid);
       end
      
       x        = x0:dx:xn;   
       lev      = fscanf(fid,'%f',[ncbot inf]);
       lev      = lev';
       fclose(fid)
      
       [bottom] = [x',bot,lev,level0.*ones(size(x'))];   

   else
       %water level
       mwl      = 0;
       x        = x0:dx:xn;
       swl      = mwl.*ones(size(x'));
       [bottom] = [x',bot,swl,level0.*ones(size(x'))];   
   end
end

%output
%=========================================================================
observ.data = obs;        %observations
observ.bot  = bottom;     %bottom file


function [model] = a11rdmod(k,fname,log,imod)

%reads data from wave model in
%.sp1 and .tab file

%extention of the wave model output file
ext      = {'.tab','.sp1'};

%number of tests included in case
ncases   = 1;

%number of columns in wave model output files
num_col  = [9,14];

%total number of headerlines in wave model output files
num_lin  = [7,7];

%check platform UNIX or DOS    
if(k == 1)
    str = ['cd ../model_io'];
else
    str = ['cd ..\model_io'];
end

eval(str);

fname = char(fname(1:6));

for jj = 1:length(ext)
    for ii = 1:ncases
       num_cas = int2str(ii);    %numbering of the case
       num_mod = int2str(imod);  %numbering of the model
   
        if (length(num_cas) < 2);
            num_cas = strcat('0',num_cas);
        end   
        f1name    = strcat(fname,num_cas,num_mod,ext{jj}); 
        
        % check if wave model output file exists
        a = exist(f1name);
        if(a == 2)      
           fid = fopen(f1name,'r');
           for kk = 1:num_lin(jj)
              line = fgetl(fid);
           end
           
           modelio      = fscanf(fid,'%f',[num_col(jj) inf]);
           fclose(fid);
        else
           fprintf(log,'%s %s\n',f1name,' does not exist');
           modelio(1:num_col(jj),:) = NaN;
        end        
        output{ii,jj} = modelio';
        tab{ii}       = output{ii,1};        
        if (jj > 1)
           sp1{ii} = output{ii,2};
        end         
    end
end

header  = {'Dist [m]','Depth [m]','Hm0 [m]','Tpeak [s]',...
      'Tm01 [s]','Tm02 [s]','Fspr []','Dir [degr]','Setup [m]'};
%output
%=========================================================================
model.tab       = tab;       %data .tab file  wave model
model.sp1       = sp1;       %data .sp1 file  wave model
model.header    = header;    %data for y-labels 

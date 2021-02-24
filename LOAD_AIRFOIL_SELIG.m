
function [dataX,dataY] = LOAD_AIRFOIL_SELIG(fileName)

hdrlns = 1;
if (strcmp(fileName,'nasasc2-0714'))
    hdrlns = 3;
elseif (strcmp(fileName,'s1020'))
    hdrlns = 2;
end

fidAirfoil = fopen(['./Airfoil_DAT_Selig/' fileName '.dat']);               
dataBuffer = textscan(fidAirfoil,'%f %f','CollectOutput',1,...              
                                         'HeaderLines',hdrlns,...
                                         'Delimiter','');
dataX = dataBuffer{1}(:,1);                                                 
dataY = dataBuffer{1}(:,2);                                                 
fclose(fidAirfoil); 

% Delete any duplicate (0,0) lines (only need one)
dataArr     = [dataX dataY];
[~,ia,~]    = unique(dataArr,'rows','stable');                              
i           = true(size(dataArr,1),1);                                      
i(ia)       = false;                                                        
dataArr(i,:)= [];                                                           
dataX       = dataArr(:,1);                                               
dataY       = dataArr(:,2);                                               

% Close the airfoil by adding final point
if (dataY(1) ~= dataY(end))                                                 
    dataX(end+1) = dataX(1);                                                
    dataY(end+1) = dataY(1);                                                
end

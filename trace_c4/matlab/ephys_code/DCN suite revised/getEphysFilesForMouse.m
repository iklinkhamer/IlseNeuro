% Load 
%
% TODO: Julius: This type of generic file loading can be handled by project-agnostic
%               code in a central library.
function fileLists =  getEphysFilesForMouse(mname)
    arguments
        mname(1,1) string
%         withBehavior(1,1) logical = false;
    end



%% --------------------------------------------------------------------
%                                   File input / Output
% ------------------------------------------------------------------------------

  disp("----------- Mouse Name ----------")
  prompt = " Please enter mouse name: e.g. ReserveMouse1 : ";
  
  if nargin > 0
      fprintf("%s%s\t(obtained from function input.)\n\n\n", prompt, mname)
  else
      mname = input(prompt,"s");
  end
  
  if isempty(mname)
      mname = "GeneralMouse";
  end
       
%   IkUtils.printIndented(["Searching for mouse  " mname]);
  
%   [res] = MouseList(mname);
%   
%   assert ...
%       ( not(res == 0) ...
%       , "Mouse not defined in pre-registered mice list." ...
%       )
%   
%   reserveMouse = res > 100;
  
%   ephysRoot = Env.getEphysRoot(reserve = reserveMouse);
  analyzedEphysDir = fullfile(ephysRoot, mname, "AnalyzedEphys");
  behavior4EphysDir = fullfile(ephysRoot, mname);
  rawEphysDir = fullfile(ephysRoot, mname);
 
  % Analyze folder contents of input and output folders
  IkUtils.printIndented("Analyzing....")
  
  filePrefix = "EphysAnalyzed";
  
  [fileNames,dateStrings,analyzedEphysFiles] = IkUtils.io.listTimestampedFiles ...
      ( filePrefix ...
      , onlyLates = false ...
      , folder = analyzedEphysDir ...
      );
  
  behaviorFileTemplateFn = @(timeStamp) ...
      fullfile ...
          ( behavior4EphysDir ...
          , sprintf("%s_%d", mname, timeStamp) ...
          , sprintf("EphysSession_%s_%d.mat", mname, timeStamp) ...
          );
  
%       behaviorFiles = arrayfun ...
%             ( @(timeStamp) dir(behaviorFileTemplateFn(timeStamp)) ...
%             , timeStamps ...
%             );
  
  fileLists = struct ...
      ( analyzedEphys = analyzedEphysFiles ...
      , rawEphys = rawEphysDir ...
      );
%   , behavior = behaviorFiles ...
    
end


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
function [DATA] = XYZProj2Mat(direction)
Name = ["X","Y","Z"];

i = direction;
Target = strcat(Name(i),"Proj.csv");
ID = fopen(Target,"r");
 if(ID == -1)
  error("file doesn't exist.")
 else 
      tline = fgetl(ID);
      c = strsplit(string(tline),',');
      while(strcmp(c(1),"Area")==0)    
         tline = fgetl(ID); 
         c = strsplit(string(tline),','); %% string(char)でcharをStringに変換
      end
   DATA = str2double(c(2));
 end
 fclose(ID);
end
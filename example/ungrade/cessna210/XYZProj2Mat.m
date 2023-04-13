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
      while(strcmp(c(1),"Comp_Areas")==0)    
         tline = fgetl(ID); 
         c = strsplit(string(tline),','); %% string(char)でcharをStringに変換
      end
   DATA = str2double(c(2:end));
 end
 fclose(ID);
end
function [PC] = PANAIRCp(proid)
wd = pwd;
pandir = strcat(wd,'\PANAIR',num2str(proid));
cd(pandir)
%panair.outから
fp = fopen('panair.out');

%データをまとめたもの
	while(1)
		buff = fgetl(fp);
        if(strncmp(buff,'                                                            panel data',70))
            break;
        end
        if buff == -1
            error('PANAIR IS ABORTED')
        end
    end
    for i = 1:6
        buff = fgetl(fp);
    end
    
    i=1;
    while(1)
        buff = fscanf(fp,'%s',3);
        DATA.PN(i,1) = fscanf(fp,'%d',1);%panel no
        buff = fscanf(fp,'%s',3);
        DATA.NN(i,1) = fscanf(fp,'%d',1);%Network no
        buff = fscanf(fp,'%s',2);
        DATA.row(i,1) = fscanf(fp,'%d',1);%row no
        buff = fscanf(fp,'%s',2);
        DATA.col(i,1) = fscanf(fp,'%d',1);%col no
        for iter = 1:8
            buff = fgetl(fp);
        end
        buff = fscanf(fp,'%s',1);
        DATA.r0(i,:) = fscanf(fp,'%f',3);%r0
        buff = fgetl(fp);buff = fscanf(fp,'%s',1);
        DATA.n0(i,:) = fscanf(fp,'%f',3);%r0
        for iter = 1:7
            buff = fscanf(fp,'%s',1);
            buff = fscanf(fp,'%f',1);
        end
        buff = fscanf(fp,'%s',1);
        DATA.area(i,1) = fscanf(fp,'%f',1);%area
        buff = fgetl(fp);
        buff = fgetl(fp);
        if(strncmp(buff,' nsngn',6))
            break;
        end
        buff = fgetl(fp);buff = fgetl(fp);buff = fgetl(fp);
        i = i+1;
        if buff == -1
            error('PANAIR IS ABORTED')
        end
    end
    while(1)
		buff = fgetl(fp);
        if(strncmp(buff,'0*b*solution',12))
            break;
        end
        if buff == -1
            error('PANAIR IS ABORTED')
        end
    end
    while(1)
        buff = fgetl(fp);
        if(strncmp(buff,'1',1))
            break;
        end
        if buff == -1
            error('PANAIR IS ABORTED')
        end
    end
    for iter = 1:7
        buff = fgetl(fp);
    end
    %ルーチン
    j=1;
    i=1;
    while(1)
        databuff = fscanf(fp,'%f',49);buff = fgetl(fp);
        DATA.Cp2ndd{j}(i,1) = databuff(48);
        buff = fgetl(fp);
        if(strncmp(buff,'1',1))
              buff = fgetl(fp);
              if(strncmp(buff,'0*b*for-mom',11))
                  while(1)
                      buff = fgetl(fp);
                      if(strncmp(buff,'0*e*for-mom',11))
                          break;
                      end
                      if buff == -1
                        error('PANAIR IS ABORTED')
                      end
                  end
                  buff = fgetl(fp);buff = fgetl(fp);
                  if(strncmp(buff,' **********',11))
                          break;
                  elseif(strncmp(buff,'0*b*solution',12))
                      j=j+1;
                      i=0;
                      for iter = 1:9
                          buff = fgetl(fp);
                      end
                  end
                  for iter = 1:6
                    buff = fgetl(fp);
                  end
              else
                for iter = 1:6
                    buff = fgetl(fp);
                end
              end
        end
        i=i+1;
        if buff == -1
            error('PANAIR IS ABORTED')
        end
    end
    fclose(fp);
    %データ整形
    for i = 1:size(DATA.Cp2ndd{1}(:,1),1)
        for jter = 1:j
            PC{DATA.NN(i,1)}.Cp{jter}(DATA.row(i,1),DATA.col(i,1)) = DATA.Cp2ndd{jter}(i,1);
        end
        PC{DATA.NN(i,1)}.X(DATA.row(i,1),DATA.col(i,1)) = DATA.r0(i,1);
        PC{DATA.NN(i,1)}.Y(DATA.row(i,1),DATA.col(i,1)) = DATA.r0(i,2);
        PC{DATA.NN(i,1)}.Z(DATA.row(i,1),DATA.col(i,1)) = DATA.r0(i,3);
        PC{DATA.NN(i,1)}.nx(DATA.row(i,1),DATA.col(i,1)) = DATA.n0(i,1);
        PC{DATA.NN(i,1)}.ny(DATA.row(i,1),DATA.col(i,1)) = DATA.n0(i,2);
        PC{DATA.NN(i,1)}.nz(DATA.row(i,1),DATA.col(i,1)) = DATA.n0(i,3);
        PC{DATA.NN(i,1)}.area(DATA.row(i,1),DATA.col(i,1)) = DATA.area(i,1);
    end

    cmd = sprintf('del ffm >nul 2>&1 && del ffmf >nul 2>&1 && del fmgp >nul 2>&1 && del fort.88 >nul 2>&1 && del fslout >nul 2>&1 && del ft?? >nul 2>&1 && del iflggp >nul 2>&1 && del ispggp >nul 2>&1 && del news >nul 2>&1 && del nli??? >nul 2>&1 && del rwms?? >nul 2>&1 && exit');
    [status,result] = system(cmd);
cd(wd)


end
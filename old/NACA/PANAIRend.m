function [status result] = PANAIRend(pid,WGS,MACH,ALPHA,REFS,n_wake,n_del,now_work)
%status = 1 :実行中
%status = 0 :計算終了
%status = -1 :計算異常終了

ALPC    = 0;
BETA    =  zeros(size(ALPHA(:)'));
BETC    =  0;


SPAN    = REFS(2,2);
SREF     = REFS(2,1);
CBAR    = REFS(2,3);

%Xrefは空力平均翼弦の前縁とする
XREF    = REFS(1,1);
YREF    =REFS(1,2);
ZREF    =REFS(1,3);


pandir = strcat(now_work,'\PANAIR',num2str(pid));
cd(pandir)
%結果ファイルの読み込み
fp = fopen('panair.out','r');
i = 1;
result.info = 1;
while(1)
    buff=fgetl(fp);
    i = i+1;
    if  feof(fp)
        result.info = 0;
        break;
    end
    if strncmp(buff,' sol-no',7) || i>=10000000
        break;
    end
end
i =0;
while(1)
    buff=fgetl(fp);
    i=i+1;
    if  feof(fp)
        result.info = 0;
        break;
    end
    if strncmp(buff,' sol-no',7) || i>=10000000
        break;
    end
end
for i= 1:3
    buff=fgetl(fp);
end

if result.info == 1
    AirData=(fscanf(fp,'%f',[13,min([size(ALPHA(:)',2),size(BETA(:)',2)])]))';
    fclose(fp);
    result.mach = MACH;
    result.alpha = AirData(:,2);
    result.beta = AirData(:,3);
    result.CL = AirData(:,4);
    result.CD = AirData(:,5);
    result.CY = AirData(:,6);
    result.CMX = AirData(:,10);
    result.CMY = AirData(:,11);
    result.CMZ = AirData(:,12);
    status = 0;
else
    fclose(fp);
    disp('PANAIR IS STOPPED\n')
    result.mach = MACH;
    result.alpha = 0;
    result.beta = 0;
    result.CL = 0;
    result.CD = 100000000;
    result.CY = 0;
    result.CMX =0;
    result.CMY = 0;
    result.CMZ = 0;
    status = -1;
end
cd(now_work)
end
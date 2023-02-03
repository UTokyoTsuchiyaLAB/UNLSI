function [status now_work] = PANAIRstart(proid,WGS,MACH,ALPHA,REFS,n_wake,n_del,isSymm)
%REFS = [XREF YREF ZREF 0;
%	     SREF BREF CREF DREF
%n_wake  wakeの数
%n_del ノズルやベース面など削除するもの

%PANIN設定
ALPC    =  mean(ALPHA);
BETA    =  zeros(size(ALPHA(:)'));
BETC    =  0;


SPAN    = REFS(2,2);
SREF     = REFS(2,1);
CBAR    = REFS(2,3);

%Xrefは空力平均翼弦の前縁とする
XREF    = REFS(1,1);
YREF    =REFS(1,2);
ZREF    =REFS(1,3);

%wgsを開いてnetnumを境界条件に
j=1;
fp = fopen(WGS,'r');
while(1)
    %'まで捨てる
    while(1)
        buff = fgetl(fp);
        if buff == -1 
            break;
        end
        if strcmp(buff(1),'''') == 1
            break;
        end
    end
    if buff == -1
        break;
    end
    BOUN(j) = fscanf(fp,'%d',1);
    j=j+1;
end
fclose(fp);

%PANIN & PANAIR実行
now_work = pwd;

%------------------------------------------
%並列実行の準備
%------------------------------------------
NUM_MACH = size(MACH(:),1);


        pandir = strcat(now_work,'\PANAIR',num2str(proid));
        cd(pandir)

        AUX = 'pin.aux'; 
        fp = fopen(AUX,'wt');
        fprintf(fp,'WGS %s\n',WGS);
        fprintf(fp,'MACH %f\n',MACH);
        fprintf(fp,'ALPHA ');
        for i = 1:size(ALPHA(:)',2)
            fprintf(fp,'%f ',ALPHA(i));
        end
        fprintf(fp,'\n');
        fprintf(fp,'ALPC %f\n',ALPC);
        fprintf(fp,'BETA ');
        for i = 1:size(BETA(:)',2)
            fprintf(fp,'%f ',BETA(i));
        end
        fprintf(fp,'\n');
        fprintf(fp,'BETC %f\n',BETC);
        fprintf(fp,'SPAN %f\n',SPAN);
        fprintf(fp,'CBAR %f\n',CBAR);
        fprintf(fp,'SREF %f\n',SREF);
        fprintf(fp,'XREF %f\n',XREF);
        fprintf(fp,'YREF %f\n',YREF);
        fprintf(fp,'ZREF %f\n',ZREF);
        fprintf(fp,'IGEOMP 1\n');
        fprintf(fp,'SYMM %d\n',isSymm);
        fprintf(fp,'EAT 0.001\n');
        fprintf(fp,'BOUN\n');
        for i=1:size(BOUN,2)
            fprintf(fp,'%d ',BOUN(i));
        end
        fprintf(fp,'\n');
        fclose(fp);

        %panin.exeと接続
        clear cmd
        cmd = sprintf('del ffm >nul 2>&1 && del ffmf >nul 2>&1 && del fmgp >nul 2>&1 && del fort.88 >nul 2>&1 && del fslout >nul 2>&1 && del ft?? >nul 2>&1 && del iflggp >nul 2>&1 && del ispggp >nul 2>&1 && del news >nul 2>&1 && del nli??? >nul 2>&1 && del rwms?? >nul 2>&1&& exit');
        [status,result] = system(cmd);

        fp = fopen('panin.inp','wt');
        fprintf(fp,AUX);
        fclose(fp);
        wd = pwd;
        clear cmd
        cmd = sprintf('cd %s && copy %s %s && panin.exe < panin.inp > panin.out && echo execute > endflag.txt && exit',wd,strcat(now_work,'\\',WGS),strcat(wd,'\\'));
        [status,result] = system(cmd);

        %a502.inに$forcesを追加
        fp = fopen('a502.in','rt');
        i=1;
        while(1)
            a502buff{i}=fgets(fp);
            i = i+1;
            if a502buff{i-1}==-1
                break;
            end
        end
        fclose(fp);

        fp = fopen('a502.in','wt');
        i=1;
        while(1)
            fprintf(fp,'%s',a502buff{i});
            i=i+1;
            if i==size(a502buff,2)-1
                break;
            end
        end
        fprintf(fp,'$forces and moments summary for nonwake networks\n');
        fprintf(fp,'*network selection for summary\n');
        fprintf(fp,'     %.1f     %.1f\n',size(BOUN,2)-n_wake-n_del,1.0);
        for i=1:size(BOUN,2)-n_wake-n_del
            fprintf(fp,'      %.1f\n',i);
        end
        fprintf(fp,'$end\n');
        fclose(fp);
        cd(now_work)


        pandir = strcat(now_work,'\PANAIR',num2str(proid));
        cd(pandir)
        clear cmd
        cmd = sprintf('del panair.out >nul 2>&1 && del ffm >nul 2>&1 && del ffmf >nul 2>&1 && del fmgp >nul 2>&1 && del fort.88 >nul 2>&1 && del fslout >nul 2>&1 && del ft?? >nul 2>&1 && del iflggp >nul 2>&1 && del ispggp >nul 2>&1 && del news >nul 2>&1 && del nli??? >nul 2>&1 && del rwms?? >nul 2>&1 && exit');
        [status,result] = system(cmd);
        fp = fopen('panair.inp','wt');
        fprintf(fp,'a502.in\n');
        fclose(fp);
        [status,result] = system('panair.exe < panair.inp');
        cd(now_work)

status = 1;
function [result PC PE] = hyper(wgsname,MACH,ALPHA,BETA,REFS,cmethods,emethods,isSym)

%wd = pwd;
%hyperdir = strcat(wd,'\HYPER',num2str(n_process));
%clear cmd
%cmd = sprintf('cd %s && copy %s %s && exit',wd,strcat(wd,'\\',wgsname),strcat(hyperdir,'\\',wgsname));
%[status,result] = system(cmd);

%cd(hyperdir)

%hyperと接続して超音速の空力性能を調べる
%REFS = [XREF YREF ZREF 0 ;
%	     SREF BREF CREF DREF

 %~ COMPRESSION METHODS
    %~ 1  Modified Newtonian
    %~ 2  Newtonian-Prandtl-Meyer
    %~ 3  Tangent Wedge
    %~ 4  Tangent Wedge Infinite Mach
    %~ 5  Old Tangent Cone
    %~ 6  Cone At Angle Of Attack (later)
    %~ 7  VanDyke Unified
    %~ 8  Blunt Body Viscous (later)
    %~ 9  Shock Expansion (later)
   %~ 10  Free Molecular Flow (later)
   %~ 11  Input value of CpStag
   %~ 12  Hankey Flat Surface
   %~ 13  Smyth Delta Wing
   %~ 14  Modified Dahlem-Buck
   %~ 15  BlastWave (later)
   %~ 16  OSUBluntBody
   %~ 17  Tangent Cone (Edwards)

 %~ EXPANSION METHODS
    %~ 1  Cp=0
    %~ 2  NewtonianPrandtlMeyer
    %~ 3  PrandtlMeyer
    %~ 4  ConeAtAngleOfAttack(later)
    %~ 5  VanDykeUnified
    %~ 6  Vacuum
    %~ 7  Shock Expansion (later)
    %~ 8  Input Value
    %~ 9  Free Molecular Flow(later)
   %~ 10  Modified Dahlem-Buck
   %~ 11  ACMempirical(later)
   %~ 12  half Prandtl-Meyer from freestream

%wgsnameよりファイル名を取得
[wgspath name ext]=fileparts(wgsname);

%CpStagの入力
CpStag = 0;

%インプットファイルを作成
fp=fopen(strcat(name,'.inp'),'wt');
fprintf(fp,'&hyp  title="made by hyper.m",\n');
fprintf(fp,'  wgsFileName="%s",\n',strcat(name,'.wgs'));
fprintf(fp,'  cmethods=');
for i = 1:size(cmethods(:),1)
    fprintf(fp,'%d,',cmethods(i));
end
fprintf(fp,'\n');
fprintf(fp,'  emethods=');
for i = 1:size(emethods(:),1)
    fprintf(fp,'%d,',emethods(i));
end
fprintf(fp,'\n');
fprintf(fp,'  mach=%.5f, sref=%.5f, cbar=%.5f,\n',MACH,REFS(2,1),REFS(2,3));
fprintf(fp,'  alpha=');
for i = 1:size(ALPHA(:),1)
    fprintf(fp,'%.5f,',ALPHA(i));
end
fprintf(fp,'  beta=');
for i = 1:size(BETA(:),1)
    fprintf(fp,'%.5f,',BETA(i));
end
fprintf(fp,'\n');
fprintf(fp,'  xref=%.5f, span=%.5f,\n',REFS(1,1),REFS(2,2));
fprintf(fp,'  yref=%.5f, zref=%.5f\n',REFS(1,2),REFS(1,3));
fprintf(fp,'  cpStag=%.5f/\n',CpStag);
fclose(fp);

%hyper.exeと接続
wd = pwd;
clear cmd
fp = fopen('hyper.min','wt');
fprintf(fp,strcat(name,'.inp\n'));
fclose(fp);
cmd = sprintf('cd %s  && hyper.exe < hyper.min > hyper.mout && exit',wd);
[status,result] = system(cmd);
if status ~=0
    sprintf('hyper is termnated abnormally');
end

fp = fopen('hyper.out','r');
i = 1;
result.info = 1;
for n=1:2
    while(1)
        buff=fgetl(fp);
        i=i+1;
        if  feof(fp)
            result.info = 0;
            break;
        end
        if strncmp(buff,' SOLUTIONS',10) == 1
            break
        end
    end
    i = 0;
end
buff=fgetl(fp);
AirData=(fscanf(fp,'%f',[6,max([size(ALPHA(:),1),size(BETA(:),1)])]))';
result.mach = MACH;
result.alpha = AirData(:,2);
result.beta = AirData(:,3);
if isSym == 1
    result.CL = AirData(:,4).*2;
    result.CD = AirData(:,5).*2;
    result.CMY = AirData(:,6).*2;
else
    result.CL = AirData(:,4);
    result.CD = AirData(:,5);
    result.CMY = AirData(:,6);
end
result.info = 1;
fclose(fp);
%{
delete('hyper.out');
delete('hyper.mout');
delete('hyper.min');
%}
%cd ..

if nargout >=2
    %hyperの圧力分布の結果を読み込んで出力する
    filename = 'hyper.dbg';

    %ファイルオープン
    fp = fopen(filename,'r');
    
    
    %要らない行を削除
    while(1)
        buff = fgetl(fp);
        if strncmp(buff,' NETWORK',8)
            break;
        end
    end

    netnum = 1;
    buff = fgetl(fp);
    while(1)
        while(1)
           databuff = fscanf(fp,'%f\n',[1,5]);
           if isempty(databuff)
               buff = fgetl(fp);
               break;
           end
           PE{netnum}.X(databuff(1),databuff(2)) = databuff(3);
           PE{netnum}.Y(databuff(1),databuff(2)) = databuff(4);
           PE{netnum}.Z(databuff(1),databuff(2)) = databuff(5);
        end
        PE{netnum}.nR = size(PE{netnum}.X,1);
        PE{netnum}.nC = size(PE{netnum}.X,2);
        buff = fgetl(fp);
        if strncmp(buff,'   # net row',12);
            break;
        end
        netnum = netnum + 1;
    end
    
    %PANELdata
    while(1)
       databuff = fscanf(fp,'%f\n',[1,11]);
       if isempty(databuff)
           buff = fgetl(fp);
           break;
       end
       PC{databuff(2)}.X(databuff(3),databuff(4))=databuff(5);
       PC{databuff(2)}.Y(databuff(3),databuff(4))=databuff(6);
       PC{databuff(2)}.Z(databuff(3),databuff(4))=databuff(7);
       PC{databuff(2)}.nx(databuff(3),databuff(4))=databuff(8);
       PC{databuff(2)}.ny(databuff(3),databuff(4))=databuff(9);
       PC{databuff(2)}.nz(databuff(3),databuff(4))=databuff(10);
       PC{databuff(2)}.area(databuff(3),databuff(4))=databuff(11);
    end
    %Cpdataまで捨てる
    nCp = 1;
    while(1)
        while(1)
            buff = fgetl(fp);
            if strncmp(buff,'PRESSURE COEFF AT',17);
                anlybuff = sscanf(buff,'%*s %*s %*s %*s %f %*s %*s %f',[1 2]);
                ALPHA(nCp) = anlybuff(1);
                BETA(nCp) = anlybuff(2);
                buff = fgetl(fp);
                break;
            end
            if buff == -1
                break;
            end
        end
        if buff == -1
            break;
        end

        %Cpの読込
        for i = 1:size(PE,2)
            Cp1{nCp} = zeros(PE{i}.nR,PE{i}.nC);
            Cp2{nCp} = zeros(PE{i}.nR,PE{i}.nC);
            Cp3{nCp} = zeros(PE{i}.nR,PE{i}.nC);
            Cp4{nCp} = zeros(PE{i}.nR,PE{i}.nC);
            for j = 1:PE{i}.nC-1
                for k = 1:PE{i}.nR-1
                    Cpbuff = fscanf(fp,'%f',[1 10]);
                    PC{i}.Cp{nCp}(k,j) = Cpbuff(8);
                    %NETWORK{i}.xCp{nCp}(k,j) = Cpbuff(5);
                    %NETWORK{i}.yCp{nCp}(k,j) = Cpbuff(6);
                    %NETWORK{i}.zCp{nCp}(k,j) = Cpbuff(7);
                    Cp1{nCp}(k,j) = Cpbuff(8);
                    Cp2{nCp}(k+1,j) = Cpbuff(8);
                    Cp3{nCp}(k,j+1) = Cpbuff(8);
                    Cp4{nCp}(k+1,j+1) = Cpbuff(8);
                end
            end
            PE{i}.Cp{nCp}=(Cp1{nCp}+Cp2{nCp}+Cp3{nCp}+Cp4{nCp})./4;
            PE{i}.Cp{nCp}(1,1) = PE{i}.Cp{nCp}(1,1).*4;
            PE{i}.Cp{nCp}(1,end) = PE{i}.Cp{nCp}(1,end).*4;
            PE{i}.Cp{nCp}(end,1) = PE{i}.Cp{nCp}(end,1).*4;
            PE{i}.Cp{nCp}(end,end) = PE{i}.Cp{nCp}(end,end).*4;
            PE{i}.Cp{nCp}(2:end-1,1) = PE{i}.Cp{nCp}(2:end-1,1).*2;
            PE{i}.Cp{nCp}(1,2:end-1) = PE{i}.Cp{nCp}(1,2:end-1).*2;
            PE{i}.Cp{nCp}(end,2:end-1) = PE{i}.Cp{nCp}(end,2:end-1).*2;
            PE{i}.Cp{nCp}(2:end-1,end) = PE{i}.Cp{nCp}(2:end-1,end).*2;
        end
        nCp = nCp+1; 
    end
    nCp = nCp-1;
    fclose(fp);
end
    

end
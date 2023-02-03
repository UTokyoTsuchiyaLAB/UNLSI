function info = wgs_make(wgsname,netnum,netname,add,dataX,dataY,dataZ,O,flip,HR)
	dataX = dataX.*HR-O(1);
	dataY = dataY.*HR-O(2);
	dataZ = dataZ.*HR-O(3);
	
	if flip(1) ~=0
		dataX=flipud(dataX);
		dataY=flipud(dataY);
		dataZ=flipud(dataZ);
	end
	if flip(2) ~=0
		dataX=fliplr(dataX);
		dataY=fliplr(dataY);
		dataZ=fliplr(dataZ);
	end
	if add == 1
		fp = fopen(wgsname,'wt');
		fprintf(fp,'Created by Naoto Morita from make_wgs\n');
	else
		fp = fopen(wgsname,'at+');
	end
	fprintf(fp,'\''%s\''\n',netname);
	setup_WGS = [netnum size(dataX,2) size(dataX,1) 0 0 0 0 0 0 0 1 1 1 0];
	fprintf(fp,'%d %d %d %d   %d %d %d   %d %d %d    %d %d %d  %d \n',setup_WGS(1),setup_WGS(2),setup_WGS(3),setup_WGS(4),setup_WGS(5),setup_WGS(6),setup_WGS(7),setup_WGS(8),setup_WGS(9),setup_WGS(10),setup_WGS(11),setup_WGS(12),setup_WGS(13),setup_WGS(14));
    for j_iter = 1:size(dataX,2)
        for i_iter =1:size(dataX,1)
            fprintf(fp,'%.10f %.10f %.10f \n',dataX(i_iter,j_iter),dataY(i_iter,j_iter),dataZ(i_iter,j_iter));
	    end
	end
	info = fclose(fp);
end
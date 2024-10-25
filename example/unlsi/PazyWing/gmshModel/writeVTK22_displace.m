%%%
%%% VTK :: writeVTK
%%%
function displacement = writeVTK22_displace(fid,NDIM,NN,NDOF,numnp,nume,node,element,U,e_stress)
displacement=zeros([numnp 3]);
%%% VTK_CELL_TYPE 
	vtk_cell_type=3;
%%% NVARS
	NVARS=1; 
	vname=char("S11");
%%% Header 
	fprintf(fid,"# vtk DataFile Version 2.0\n");
	fprintf('DataName\n');
	fprintf(fid,"ASCII\n");
%%% Node 
	fprintf(fid,"DATASET UNSTRUCTURED_GRID\n");
	fprintf(fid,"POINTS %d float\n",numnp);
	for i=1:numnp
		xc(1)=0.0;xc(2)=0.0;xc(3)=0.0;
		for j=1:NDIM
			xc(j)=node(i,j);
            if( NDIM == 2)
                xc(1) = xc(1)+U(2*i-1);
                xc(2) = xc(2)+U(2*i);
			    % fprintf(fid," %e %e 0.0\n",U(2*i-1),U(2*i));
            end
		end
		fprintf(fid,"%e %e %e\n",xc(1),xc(2),xc(3));
        x_displace = [xc(1) xc(2) xc(3)];
        displacement(i,:) = x_displace;
	end
%%% Element
	fprintf(fid,"CELLS %d %d\n",nume,nume*(1+NN));
	for i=1:nume
		for j=1:NN
			nid(j)=element(i,j);
		end
		fprintf(fid,"%d ",NN);
		for j=1:NN
			fprintf(fid,"%d ",nid(j)-1);
		end
		fprintf(fid,"\n");
	end
%%% type
	fprintf(fid,"CELL_TYPES %d\n",nume);
    for i=1:nume
	    fprintf(fid,"%d\n",vtk_cell_type);
    end
%%% Result for NODE
	fprintf(fid,"POINT_DATA %d\n",numnp);
%%% Deformation 
	fprintf(fid,"VECTORS DISPLACEMENT float\n");
	for i=1:numnp
		if( NDIM == 2)
			fprintf(fid," %e %e 0.0\n",U(2*i-1),U(2*i));
        end
	end
%%% Result for ELEMENT
	fprintf(fid,"CELL_DATA %d\n",nume);
%%% Element Stress 
	for i=1:NVARS
		fprintf(fid,"SCALARS %s float 1\n",vname(i,:));
		fprintf(fid,"LOOKUP_TABLE default\n");
		for j=1:nume
			fprintf(fid,"%e\n",e_stress(j));
		end
	end
end
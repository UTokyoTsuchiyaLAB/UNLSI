clear;
%%%
%%% NDIM definitions: Dimension of problem
%%%
NDIM=2;%要素の次元
%%%
%%% NDOF definitions: Number of DOFs per node
%%%
NDOF=2;%接点当たりの自由度
%%%
%%% NN definitions: Number of nodes per element
%%%
NN=2;%接点数
% figure;
N = 10;
airfoil = CST_airfoil([-0.1294 -0.0036 -0.0666], [0.206 0.2728 0.2292],0,N);

lowerSurf = flip(airfoil(1:N/2,:));
upperSurf = airfoil(N/2+2:N+1,:);
leadingEdge = airfoil(N/2+1,:);
trailingEdge = airfoil(1,:);
% disp(lowerSurf);
% disp(upperSurf);

%%%%接点の指定
numnp = N+(N-4)/2;
% numnp = N;
node = zeros([numnp 2]);
node(1,:)=leadingEdge;
for i =1:(N-2)/2
    node(2*i,:) = upperSurf(i,:);
    node(2*i+1,:) = lowerSurf(i,:);
end
node(N,:) = trailingEdge;
for j=1:(N-4)/2
    nodex = (upperSurf(j,:)+lowerSurf(j+1,:))/2;

    node(N+j,:) = nodex;
end

%%%要素の指定
nume = N + (N-2)/2+ (N-4)*2; 
element = zeros([nume 2]);
element(1,:)=[1 2];
for k=1:((N-2)/2)
    element(2*k,:)=[2*k 2*k+2];%上面
    element(2*k+1,:)=[2*k-1 2*k+1];%下面
    element(N+k,:)=[2*k 2*k+1];%縦
end
element(N,:)=[N-1 N];
for  l =1:(N-4)/2
    element(N+k+4*l-3,:)=[2*l N+l];
    element(N+k+4*l-2,:)=[2*l+1 N+l];
    element(N+k+4*l-1,:)=[2*l+2 N+l];
    element(N+k+4*l,:)=[2*l+3 N+l];
end

scatter(node(:,1),node(:,2),'cyan');
hold on
for m=1:length(element(:,1))
    x=[node(element(m,1),1) node(element(m,2),1)];
    y=[node(element(m,1),2) node(element(m,2),2)];
    plot(x,y,'blue');
end

xlim([0,1]);
ylim([-0.5,0.5]);
E=100;%100.0e3;
A=10000;%10000.0;
nspc=3;
spc(1,:)=[11 1];
spc(2,:)=[11 2];
% spc(1,:)=[116 1];
% spc(2,:)=[116 2];
spc(3,:)=[1 1];
% spc(4,:)=[1 2];

ncload=(N-2)/2;
cload = zeros([m 2]);
vcload = zeros([m 1]);
for m=1:ncload
    cload(m,:)=[2*m 2];
    vcload(m)= -100;
end

ndload = 0;

%%%
%%% Global stiffness (NDOF*numnp) x (NDOF*numnp)
%%%
KG=zeros(NDOF*numnp,NDOF*numnp);
%%%
%%% Loop over elements
%%%
for ie=1:nume
%%% Local element coordinate
	for i=1:NN
		for j=1:NDIM
			nodeL(i,j)=node(element(ie,i),j);
		end
	end
%%% Local elememt stiffness element
	kL=k_truss2d2(nodeL,E,A);
%%% Add element 1 to Global stiffness 
	for i=1:NN
		for ii=1:NDOF
			irow=NDOF*(element(ie,i)-1)+ii;
			for j=1:NN
				for jj=1:NDOF
					icol=NDOF*(element(ie,j)-1)+jj;
					KG(irow,icol)=KG(irow,icol)+kL(NDOF*(i-1)+ii,NDOF*(j-1)+jj);
				end
			end
		end
	end
end
%%%
%%% Load vector (condentrated load)
%%%
FC=zeros(NDOF*numnp,1);
for i=1:ncload
	npos=NDOF*(cload(i,1)-1)+cload(i,2)
	FC(npos)=FC(npos)+vcload(i);
end
%%%
%%% Load vector (distributed load)
%%%
FD=zeros(NDOF*numnp,1);
for i=1:ndload
	ie=dload(i);
	for j=1:NN
		for k=1:NDIM
			nodeL(j,k)=node(element(ie,j),k);
		end
	end
	bf=bf_truss2d2(nodeL,A,vdload(i,:));
	for i=1:NN
		for ii=1:NDOF
			irow=NDOF*(element(ie,i)-1)+ii;
			FD(irow)=FD(irow)+bf(NDOF*(i-1)+ii);
		end
	end
end
%%%
%%% TOTAL LOAD
%%%
	F=FC+FD;
%%%
%%% constrained conditions
%%%
count=1;
for i=1:nspc
	npos=NDOF*(spc(i,1)-1)+spc(i,2);
	row(count,:)=KG(npos,:);
	KG(npos,:)=0;
	KG(:,npos)=0;
	KG(npos,npos)=1;
	F(npos)=0.0;
	count=count+1;
end
%%%
%%% computation 
%%%
fprintf("### Computation started\n");
U=KG\F;
fprintf("### Computation finished\n");
fprintf("### Result displacement\n")
for i=1:numnp
	fprintf("%d ",i);
	for j=1:NDOF
		fprintf("%e ",U(NDOF*(i-1)+j));
	end
	fprintf("\n");
end
%%% post processing :: stress
%%% for element stress/S11
e_stress=zeros(nume);
fprintf("Element force\n")
for ie=1:nume
	for j=1:NN
		for k=1:NDIM
			nodeL(j,k)=node(element(ie,j),k);
		end
		for k=1:NDOF
			dispL(j,k)=U(NDOF*(element(ie,j)-1)+k);
		end
	end
	s11=s_truss2d2(nodeL,E,A,dispL);
	fprintf("%d %e\n",ie,s11);
	e_stress(ie)=s11;
end
%%% post processing :: reaction force
fprintf("Reaction force\n")
count=1;
for i=1:nspc
	f=-FD(NDOF*(spc(i,1)-1)+spc(i,2));
	for k=1:NDOF*numnp
		f=f+row(count,k)*U(k);
	end
	count=count+1;
	fprintf("%d %d %e\n",spc(i,1),spc(i,2),f);
end
%%%
%%% for visualization
%%%
	fid=fopen('naca4412.vtk','w');
	displaced = writeVTK22_displace(fid,NDIM,NN,NDOF,numnp,nume,node,element,U,e_stress);
	fclose(fid);
%%%
%%% finished
%%%
scatter(displaced(:,1),displaced(:,2))
for m=1:length(element(:,1))
    x=[displaced(element(m,1),1) displaced(element(m,2),1)];
    y=[displaced(element(m,1),2) displaced(element(m,2),2)];
    plot(x,y,'magenta');
end
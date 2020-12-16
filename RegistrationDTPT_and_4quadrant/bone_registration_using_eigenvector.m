
clear all 

ignore_proximal=false;
ignore_distal=false;

%if any one proximal or distal part is missing than change false to true. 


if ignore_distal==true
     DT=[];
else
    distal_data_path='Nuclei_and_Cells_DT_S84_m1_mut/';    
    [numbers,txt,raw] = xlsread([distal_data_path,'Tile_coordinates.xlsx']);
       coordinates = zeros(size(txt,1)-3,5);
        for i = 4:size(txt,1),
            temp =  char(txt(i,1));
            res = strsplit(temp,'_POS');
            coordinates(i-3,1) = str2num(char(res(2)));
            coordinates(i-3,2:5) = numbers(i-3,:);
        end
        tile=coordinates(:,2:end);
   
    
    %Read (Characteristics) files
    start=0;
     for position = coordinates(:,1)'      
        a=load(strcat(distal_data_path,'c_n_pos',num2str(position),' (Characteristics).mat'));
        fi = find(coordinates(:,1) == position);
        finish=size(a.G.nuc.centroids,1);
        DT(start+1:start+finish,:)=a.G.nuc.centroids+repmat([tile(fi,2), -tile(fi,1),tile(fi,3)],finish,1);
        start=start+finish;
    end
    
   
    
end



if ignore_proximal==true
      PT=[];
else
      proximal_data_path='Nuclei_and_Cells_PT_S84_m1_mut/';
      [numbers,txt,raw] = xlsread([proximal_data_path,'Tile_coordinates.xlsx']);
       coordinates = zeros(size(txt,1)-3,5);
        for i = 4:size(txt,1),
            temp =  char(txt(i,1));
            res = strsplit(temp,'_POS');
            coordinates(i-3,1) = str2num(char(res(2)));
            coordinates(i-3,2:5) = numbers(i-3,:);
        end
        tile=coordinates(:,2:end);
   

        %Read (Characteristics) files
        start=0;
        for position = coordinates(:,1)'      
            a=load(strcat(proximal_data_path,'c_n_pos',num2str(position),' (Characteristics).mat'));
            fi = find(coordinates(:,1) == position);
            finish=size(a.G.nuc.centroids,1);
            PT(start+1:start+finish,:)=a.G.nuc.centroids+repmat([tile(fi,2), -tile(fi,1),tile(fi,3)],finish,1);
            start=start+finish;
        end
        
end


indDT=1:size(DT,1);
DTPT=[DT;PT];
indPT=size(DT,1)+1:size(DTPT,1);
DTPT=DTPT-mean(DTPT);

h1=figure;
shp = alphaShape(DTPT);
plot(shp)
[tetrahedron,V] = alphaTriangulation(shp);
trep = triangulation(tetrahedron, V);
[tri, Xb] = freeBoundary(trep);
V=Xb;
[vec,val]=eig(cov(V));
newDTPT=DTPT*vec;
alphaShapeAll=V;


%DT
if ignore_distal==true
    alphaShapeDT=[];
else
	cDT=DT-mean(DT);
	shp = alphaShape(cDT);
	[tetrahedron,V] = alphaTriangulation(shp);
	trep = triangulation(tetrahedron, V);
	[tri, Xb] = freeBoundary(trep);
	alphaShapeDT=Xb;
end

%PT
if ignore_proximal==true
     alphaShapePT=[];
else
	cPT=PT-mean(PT);
	shp = alphaShape(cPT);
	[tetrahedron,V] = alphaTriangulation(shp);
	trep = triangulation(tetrahedron, V);
	[tri, Xb] = freeBoundary(trep);
	alphaShapePT=Xb;
end


  

h2=figure; 
set(gcf, 'PaperSize', [13 16]);
set(gcf, 'PaperPosition', [0 0 13 16]);

XL=0.09;XR=0.01;XGap=0.1;Row=4;
YT=0.08;YB=0.08;YGap=0.07;Col=2;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 


for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        marray=[XPos,YPos,Width,Height];
        subplot('Position',marray);
         
        if j==1
        plot3(DTPT(indDT,1),DTPT(indDT,2),DTPT(indDT,3),'b.')
        hold on 
        plot3(DTPT(indPT,1),DTPT(indPT,2),DTPT(indPT,3),'r.')
        end
        
        if j==2
            plot3(newDTPT(indDT,1),newDTPT(indDT,2),newDTPT(indDT,3),'b.')
            hold on 
            plot3(newDTPT(indPT,1),newDTPT(indPT,2),newDTPT(indPT,3),'r.')
        end
        
        
        xlabel('x');ylabel('y');zlabel('z')
        if i==1
            view(90,90)
        elseif i==2
            view(90,0)
        elseif i==3
            view(0,0)
        end
        
        if chro==1
            title('Before registration','fontweight','normal')
        end
        if chro==2
            title('After registration','fontweight','normal')
        end
  
        XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end

    


if ignore_distal==false
   %dlmwrite(strcat(distal_data_path,'Alignment_matrix.dat'),[vec],'\t');
   saveas(h1,[distal_data_path,'BoneAlphaShape','.png']);
   saveas(h2,[distal_data_path,'registration','.png']);
   M=DTPT(indDT,:);
   save(strcat(distal_data_path,'AlignedXYZ.mat'),'M','alphaShapeAll','alphaShapePT', 'alphaShapeDT','vec','val'); 
end


if ignore_proximal==false
   %dlmwrite(strcat(proximal_data_path,'Alignment_matrix.dat'),[vec],'\t');
   saveas(h1,[proximal_data_path,'BoneAlphaShape','.png']);
   saveas(h2,[proximal_data_path,'registration','.png']);
   M=DTPT(indPT,:);
   save(strcat(proximal_data_path,'AlignedXYZ.mat'),'M','alphaShapeAll','alphaShapePT', 'alphaShapeDT','vec','val');
end


close all 

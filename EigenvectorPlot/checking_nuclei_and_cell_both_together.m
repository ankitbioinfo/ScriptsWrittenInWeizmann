

clear all 

% 'G:\Labs\EliZelzer\agrawala\Temp\SarahClusterLeftRightTibiaBone'

% opt.path={  'Nuclei_and_Cells_DT_S17_m1_mut\','Nuclei_and_Cells_DT_S18_m6_wt\',...
%             'Nuclei_and_Cells_PT_S17_m1_mut\','Nuclei_and_Cells_PT_S18_m6_wt\',...    
%             'Nuclei_and_Cells_PT_S84_m4_wt\','Nuclei_and_Cells_PT_S84_m5_mut\',...
%             'Nuclei_and_Cells_DT_S17_m2_wt\','Nuclei_and_Cells_DT_S18_m2_mut\',...
%             'Nuclei_and_Cells_DT_S84_m3_wt\','Nuclei_and_Cells_DT_S84_m1_mut\',...
%             'Nuclei_and_Cells_DT_S51_m2_wt\','Nuclei_and_Cells_DT_S84_m4_wt\',...
%             'Nuclei_and_Cells_DT_S84_m5_mut\','Nuclei_and_Cells_PT_S17_m2_wt\',...
%             'Nuclei_and_Cells_PT_S18_m2_mut\','Nuclei_and_Cells_PT_S84_m3_wt\',...
%             'Nuclei_and_Cells_PT_S84_m1_mut\','Nuclei_and_Cells_PT_S51_m2_wt\',...
%             'G:\Labs\EliZelzer\agrawala\Temp\All_wild_type\data\Nuclei_and_Cells_DU_S51_m2_wt\',...
%             'G:\Labs\EliZelzer\agrawala\Temp\All_wild_type\data\Nuclei_and_Cells_DU_S84_m2_wt\',...
%             'G:\Labs\EliZelzer\agrawala\Temp\All_wild_type\data\Nuclei_and_Cells_DU_S84_m3_wt'};

opt.path={      'Nuclei_and_Cells_DT_S17_m1_mut/','Nuclei_and_Cells_DT_S18_m6_wt/',...
            'Nuclei_and_Cells_PT_S17_m1_mut/','Nuclei_and_Cells_PT_S18_m6_wt/',...    
            'Nuclei_and_Cells_DU_S51_m2_wt/', 'Nuclei_and_Cells_DU_S84_m2_wt/', 'Nuclei_and_Cells_DU_S84_m3_wt/'};

%opt.path=  {  'Nuclei_and_Cells_DU_S84_m3_wt/'};


% 
% HypertrophicTileIndex={[19:24],[24:31],...
%                        [1:6], [1:5],...
%                        [26:31],[19:24],...
%                        [1:7], [1:5],...  
%                        [20:22],[18:20],...
%                        [19:24], [20:22],...
%                        [1:5],[1:7],...
%                        [27:29],[25:30],...
%                        [1:3],[1:8],...
%                        [1:6],[15:21],[16:21]};  
%                          



for nof=1:length(opt.path)
        path=opt.path{nof};
        [numbers,txt,raw] = xlsread([path,'Tile_coordinates.xlsx']);
        coordinates = zeros(size(txt,1)-3,5);
        for i = 4:size(txt,1),
                temp =  char(txt(i,1));
                res = strsplit(temp,'_POS');
                coordinates(i-3,1) = str2num(char(res(2)));
                coordinates(i-3,2:5) = numbers(i-3,:);
        end
        tile=coordinates(:,2:end);

        ccount=1;ncount=1;
    
        cstart=0;nstart=0;
        clear nuc_centroid 
        paulnuclei=0;
        
        index=0;        
        for position = coordinates(:,1)'
                    %[nof,position]
                    %disp(num2str(position));
                    load([path,'c_n_pos',num2str(position),' (Characteristics).mat']);
                    fi = find(coordinates(:,1) == position);                    
%                     ankit(fi,:)=[size(G.cel.centroids,1), size(G.nuc.centroids,1), G.inter.n_overlapping];
%                     ankur(fi,:)=[size(G.inter.volume_ratio), sum(G.inter.volume_ratio>1)  ];
                    indtemp = find(G.inter.volume_ratio>1);                       
%                     [~,ia]=unique(G.nuc.centroids,'rows');
%                     %paulnuclei=paulnuclei+length(ia);
%                     n1=length(ia);%n1=size(G.nuc.centroids,1);
%                     [position,size(G.nuc.centroids,1),n1]                                   
                    nfinish=length(indtemp);
                    if nfinish>0
                    %nuc_centroid(nstart+1:nstart+nfinish,:)=[index+indtemp, repmat([fi],nfinish,1)];
                    nuc_centroid(nstart+1:nstart+nfinish,:)=[G.cel.volume(indtemp), G.nuc.volume(indtemp), repmat([fi],nfinish,1)];
                    nstart=nstart+nfinish;
                    %index=index+n1;        
                    end
        end
    
        dlmwrite(strcat(path,'NucleiIndex.dat'),nuc_centroid,'\t') 
        
        a1=load([path,'all_cells_nuclei.mat']);
        b1=load([path,'all_cells.mat']);
        nuc=a1.all_cells_nuclei;
        cel=b1.all_cells;
        
        %[size(cel,1), size(nuc,1)]
end








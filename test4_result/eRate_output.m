file = 'eRate_LR_post1_1';
load([file,'.mat']); eval(['eRate = ', file,';']);
fileID = fopen([file,'.txt'],'w');

data = eRate{:,["nSub","ez_h0","rz_h0","ep_h0","rp_h0"]};
fprintf(fileID,'$1/%3u$ & %7.2E &  -- & %7.2E &  -- \\\\\n',data(1,[1,2,4]));
fprintf(fileID,'$1/%3u$ & %7.2E & %4.2f & %7.2E & %4.2f\\\\\n',data(2:end,:)');
fprintf(fileID,eRate.Properties.Description);
fclose(fileID);

%%
file = 'eRate_LR_post2_1';
load([file,'.mat']); eval(['eRate = ', file,';']);
fileID = fopen([file,'.txt'],'w');

data = eRate{:,["nSub","ep_h0","rp_h0"]};
fprintf(fileID,'$1/%3u$ & %7.2E &  -- \\\\\n',data(1,[1,2]));
fprintf(fileID,'$1/%3u$ & %7.2E & %4.2f \\\\\n',data(2:end,:)');
fprintf(fileID,eRate.Properties.Description);
fclose(fileID);
%%
file = 'eRate_LR_2';
load([file,'.mat']); eval(['eRate = ', file,';']);
fileID = fopen([file,'.txt'],'w');

data = eRate{:,["nSub","eu_h1","ru_h1","ez_h0","rz_h0","ep_h0","rp_h0"]};
fprintf(fileID,'$1/%3u$ & %7.2E &  --  & %7.2E &  --  & %7.2E &  -- \\\\\n',data(1,[1,2:2:6]));
fprintf(fileID,'$1/%3u$ & %7.2E & %4.2f & %7.2E & %4.2f & %7.2E & %4.2f\\\\\n',data(2:end,:)');
fprintf(fileID,eRate.Properties.Description);
fclose(fileID);
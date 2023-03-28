file = 'eRate_P1_1';
load([file,'.mat']); eval(['eRate = ', file,';']);
fileID = fopen([file,'.txt'],'w');

data = eRate{:,["nSub","eu_h1","ru_h1","ez_h0","rz_h0","ez_div","rz_div","ep_h0","rp_h0"]};
fprintf(fileID,'$1/%3u$ & %7.2E &  --  & %7.2E &  --  & %7.2E &  --  & %7.2E &  -- \\\\\n',data(1,[1,2:2:8]));
fprintf(fileID,'$1/%3u$ & %7.2E & %4.2f & %7.2E & %4.2f & %7.2E & %4.2f & %7.2E & %4.2f\\\\\n',data(2:end,:)');
fprintf(fileID,eRate.Properties.Description);
fclose(fileID);
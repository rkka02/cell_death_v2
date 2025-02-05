tcfFilePath = "230510.174730.HeLa_Hoechst.001.Group1.A2.T001P02.TCF";
% htData = single(permute(h5read(tcfFilePath, '/Data/3D/000000'), [2,1,3])); % orignal format
timeFrameIndex = 0;
dataPath = sprintf('/Data/3D/%06d', timeFrameIndex);

htData = ReadLDMTCFHT(tcfFilePath, timeFrameIndex);% new format
% htData = htData + (1.495-1.337);
% htData = htData .* 10000;
figure();orthosliceViewer(htData);

h5writeattr('test.TCF','/','FormatVersion','1.3','TextEncoding','UTF-8')
h5create('test.TCF','/Data/3D',[size(htData,1),size(htData,2),size(htData,3)]);
h5write('test.TCF','/Data/3D',htData);
h5writeattr('test.TCF',)


% h5disp("231107.122428.A2_1425097_CTA-GeoMX.001.Group1.A1.S001.h5")

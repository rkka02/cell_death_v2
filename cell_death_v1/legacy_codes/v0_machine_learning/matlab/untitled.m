clear all; clc; clearvars;

TCF_file=dir("*.TCF");
TCF_name=TCF_file(2).name;
TCF_info=h5info(TCF_name);
TCF_attr=h5readatt(TCF_name,'/Data/3D','resolutionx')
HT3D=h5read(TCF_name,'/Data/3D/000000');
figure, sliceViewer(HT3D)
HT2D=HT3D(:,:,30);
figure, imagesc(squeeze(HT2D)),axis image off, colormap gray
HT3D=rescale(HT3D,1.33,1.4);
HighRIgranules_threshold=multithresh(HT2D,4);
LD_mask=HT3D>1.38;
figure, sliceViewer(HT3D)
figure, sliceViewer(LD_mask)
LD_meanRI=mean(HT3D((HT3D>1.38)));
% figure, orthosliceViewer(HT2D)
figure, orthosliceViewer(HT2D>HighRIgranules_threshold(2))

Organoid_thresh=multithresh(HT3D,2);
Organoid_mask=HT3D>Organoid_thresh(2);
figure, orthosliceViewer(Organoid_mask)

morphometrics=regionprops3(Organoid_mask,'All');
figure, orthosliceViewer(morphometrics.Image{4,1})
eccentricity=morphometrics(4,:).PrincipalAxisLength(2)/morphometrics(4,:).PrincipalAxisLength(1);


function Int_ref = getmeanventricleT2(MRI_file, outfile)

MRI=SurfStatReadVol1(MRI_file);
MRIROI = MRI.data(70:350, 230:300, 130:330);
options = statset('Display','final', 'MaxIter', 100);
data = MRIROI(MRIROI>12);

obj = gmdistribution.fit(data, 3, 'Options', options);

[Mtmp bb] =sort(obj.mu);
Class=obj.cluster(data);

CSF_mean = mean(data(Class==bb(1)));
GM_mean = mean(data(Class==bb(2)));
WM_mean = mean(data(Class==bb(3)));


Class1 = MRI;
Class1.data(MRI.data<=(CSF_mean + GM_mean)/2)=1;
Class1.data(MRI.data>(CSF_mean + GM_mean)/2)=2;
Class1.data(MRI.data>(WM_mean + GM_mean)/2)=3;
Class1.file_name=strcat(outfile, '_Brainclass.mnc')
SurfStatSecureWriteVol1(Class1);
system(['gzip -f ' strcat(outfile, '_Brainclass.mnc')]);


MRIROI = MRI.data(170:280, 230:350, 150:250);
ClassROI=Class1.data(170:280, 230:350, 150:250);
Int_ref = mean(MRIROI(ClassROI==3 & MRIROI>3));
save(strcat(outfile, '_T2norm.mat'),'Int_ref')

Ventricle = MRI; 
Ventricle.data = zeros(size(MRI.data));
Ventricle.data(170:280, 230:350, 150:250) = MRI.data(170:280, 230:350, 150:250)>3;
Ventricle.data(170:280, 230:350, 150:250) = Ventricle.data(170:280, 230:350, 150:250) .* Class1.data(170:280, 230:350, 150:250)==3;
Ventricle.file_name=strcat(outfile, '_Ventricle.mnc')
SurfStatSecureWriteVol1(Ventricle);
system(['gzip -f ' strcat(outfile, '_Ventricle.mnc')]);

save(strcat(outfile, '_T2norm.mat'),'Int_ref')



%%   L-RICS function based on the following article:
%$£   
%$£  Lorenzo Scipioni, Melody Di Bona, Giuseppe Vicidomini, Alberto Diaspro & Luca Lanzanò 
%$£  "Local raster image correlation spectroscopy generates high-resolution intracellular diffusion maps" 
%$£  COMMUNICATIONS BIOLOGY 1:(2018) 1:10 
%$£
%$£
%% INPUT Menus:
%$£ 
%$£ Data           : RICS Dataset as a .tif or .avi file

%$£ Mask size      : Size of local ROI (possible values: 13,17,21,25,29,33 )
%$£ D scale        : Dmin and Dmax value in um2/s for visualization and for export in .tif
%$£ Sigma_Smooth   : Size for the Gaussian smoothing
%$£ Pixel time (us): Temporal delay between 2 consecutive pixels

%% OUTPUT:
%$£ D_Map(.tif file): diffusion Map between Dmin and Dmax scaled between 0
%                    and 256
%$£ D_Map           : diffusion Map
%$£ Phase_map_filt  : Phase Map with gaussian filter
%$£ Phase_map_raw   : Phase Map without gaussian filter

%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£ 
%$£                                                                     %$£
%$£                     Scipioni Lorenzo and Luca Lanzanò               %$£
%$£         Istituto Italiano di Tecnologia - Nanoscopy Department      %$£
%$£                      User-Friendly Version (Sep 2019)               %$£
%$£                                                                     %$£
%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£


function LRICS_simple();
close all

[filename,pathname, ~] = uigetfile({'*.tif';'*.avi'},'Please select an image stack');
filenamefull = [pathname, filename];   
if strcmp(filename(end-3:end),'.tif')
M=LRICS_readfiletif(filenamefull);  
elseif strcmp(filename(end-3:end),'.avi')
    M=AVI2Matrix(filenamefull);
end
% X=size(M,1);
% X1=64;
% sub=round((X-X1)/2);
% M=M(sub:end-sub,sub:end-sub,1:200);
T=size(M,3);
% add menu
prompt = {'Mask Size (13,17,21,25,29,33)','Pixel Time (us)','D scale in um2/s(Dmin Dmax)','Gauss Filter (px)'};
dlg_title = 'Input parameters'; 
num_lines = 1;
def = {'13','50','0 100','2'};
inputdata = inputdlg(prompt,dlg_title,num_lines,def);

MaskS=str2num(inputdata{1});
PixelTime=str2num(inputdata{2}); 
ColorScale=str2num(inputdata{3});   % repetition smooth images 
GaussAverage=str2num(inputdata{4});


path1=which('Param_Cal_6_2_16.mat');
load(path1);
mask_inp=round((MaskS-1)/2);
[val,idx]=min(abs(mask_vect-mask_inp));
mask=mask_vect(idx);
inputdata{1}=2*mask+1;
ColorScale=round(ColorScale);
TimeAverage=0;
Frame_Vect=0;
[D_map_3, ~,Other]=LRICS(M,mask,TimeAverage,GaussAverage,Frame_Vect);
Phase_map_raw=Other.PHASEx_y{1}(:,:,1);
Phase_map_filt=Other.PHASEx_y{1}(:,:,2);
CorrF=50/PixelTime;
D_map=CorrF*D_map_3{1,1};

figure
imagesc(D_map, ColorScale);
colormap hsv
colorbar
axis image
title(['D_a_v=',num2str(mean(D_map,'all'),3),'+/-',num2str(std(D_map,1,'all'),2),' um2/s']);


answer = questdlg('Save data?');
if strcmp(answer,'Yes')==1
%save data in Matlab
save([filenamefull(1:end-4),'_LRICS.mat'] ,'inputdata', 'D_map', 'Phase_map_filt', 'Phase_map_raw');
%export images (using Matlab)
A1=double(D_map);
MaxValImg=ColorScale(2);
MinValImg=ColorScale(1);
Aout=(A1-MinValImg)/(MaxValImg-MinValImg);
outputFileName = [filenamefull(1:end-4), '_Dmap_min',num2str(MinValImg,2),'max',num2str(MaxValImg),'.tif'];
delete outputFileName ;
imwrite(Aout, outputFileName);
end

end

% function read
function A=LRICS_readfiletif(fname)
% fname = [filename, '.tif'];
info = imfinfo(fname);
nslice = numel(info);
A=imread(fname, 1);  % read tif stack data
for k = 2:nslice
    B = imread(fname, k);
    A=cat(3,A,B);
end
end
function a=rgb2gray(b)
b=double(b);
b(:,:,1)=b(:,:,1)*0.299;
b(:,:,2)=b(:,:,2)*0.587;
b(:,:,3)=b(:,:,3)*0.114;
a=sum(b,3);
end
function M=AVI2Matrix(Filename)

A=VideoReader(Filename);
nFrames=A.Duration*A.FrameRate;
M=read(A,1);

if size(M,3)==3
    RGBflag=1;
    M=rgb2gray(M);
else
    RGBflag=0;
end

if RGBflag==1
for i=2:nFrames
    M(:,:,i)=rgb2gray(read(A,i));
end
else
for i=2:nFrames
    M(:,:,i)=read(A,i);
end
end    
end
function Col_bar_script(x,K,Colormap,ImgTitle,BarTitle,pos)

if nargin<4
    ImgTitle='';
    BarTitle='';
if nargin==1
    K=0;
    Colormap=jet;
end
end

if length(K)==2
m=K(1);
M=K(2);
else
M=absmax(nonzeros(x));
m=absmin(nonzeros(x));
if isempty(m)
    m=0;
end
if isempty(M)
    M=0;
end
end

if m==M
    M=m+10^-5;
end
positionVector1 = [0.05, 00.05, 0.75, 0.9];
title(BarTitle)
subplot('Position',positionVector1)
imagesc(x,[m M])
title(ImgTitle,'FontSize',6)
colormap(Colormap)
axis image
Figure_Format_script(pos)
g=gcf;
v=g.Position;

height=125/v(3)*3;
y=201;
l=flip((m:(M-m)/(y-1):M));
col_bar=repmat(l',1,2);
positionVector2 = [0.81, (1-height)/2, 0.04, height];
subplot('Position',positionVector2)
imagesc(col_bar)
%axis image
set(gca,'YAxisLocation','right');
set(gca,'YTick',(1:round(y/10):y))
set(gca,'YTickLabel',flip(round((m:(M-m)/10:M),1)))
set(gca,'XTickLabel',[])
positionVector1 = [0.05, 00.05, 0.75, 0.9];
title(BarTitle,'FontSize',15)
subplot('Position',positionVector1)
imagesc(x,[m M])
title(ImgTitle,'FontSize',11)
colormap(Colormap)
axis image
Figure_Format_script(pos)
end

function Figure_Format_script(pos)
FigHandle = gcf;
switch pos
    case 1
set(FigHandle, 'Position', [50, 250, 500, 500]);
    case 2
set(FigHandle, 'Position', [550, 250, 500, 500]);
    case 3
set(FigHandle, 'Position', [1050, 250, 500, 500]);
end
%axis image
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'Ztick',[])
set(gca,'Zticklabel',[])
xlabel('')
ylabel('')
zlabel('')
%FName=get(gcf,'FileName');
%saveas(FigHandle,FName);
%saveas(FigHandle,strcat(FName(1:end-4),'.tif'));

%close all
end


function [D_map,G0_map,Other]=LRICS(Matrix,Mask,Mask_Time,Sigma_Smooth,Frame_Vect)

switch nargin
    case 2
        Mask_Time=0;
        Sigma_Smooth=0;
        Frame_Vect=0;
    case 3
        Sigma_Smooth=0;
        Frame_Vect=0;
    case 4
        Frame_Vect=0;
end

%$£ Disciminates between carpet or stack
if size(Matrix,3)==1

%% $£ Carpet analysis
[D_map,G0_map,Int,Phase,Modulus]=L_RICS_analysis_Concentration(Matrix,Mask*2+1,Mask_Time,Frame_Vect);

%$£ Storing parameters
Other=struct;
Other.Type='Line';
Other.Intensity=Int;
Other.Phase=Phase;
Other.Modulus=Modulus;

else    
    
%% $£ Stack analysis    
D_map=cell(1);
G0_map=cell(1);
[ACFx,ACFy,gs_x,gs_y,Phase_x,Phase_y,Mod_x,Mod_y,Mod_bkg1_x,Mod_bkg1_y]=RICS_Local_ACF_MultiImages_Concentration(Matrix,Mask_Time,Mask,Sigma_Smooth,Frame_Vect);
% path=which('LRICS');
% path=path(1:end-7);
path1=which('Param_Cal_6_2_16.mat');
load(path1);

pos=find(mask_vect==Mask);

%$£ Scaling phase to diffusion
for i=1:length(Frame_Vect)
    if Sigma_Smooth==0
D_map{i}=real(-Param_Cal(pos,3)*log((Phase_x{i}(:,:,1)-Param_Cal(pos,1))/Param_Cal(pos,2)));
    else
D_map{i}=real(-Param_Cal(pos,3)*log((Phase_x{i}(:,:,2)-Param_Cal(pos,1))/Param_Cal(pos,2)));
    end
    F=(Mod_x{i}(:,:,1)./Mod_bkg1_x{i}(:,:,1))-1;
A=0.0357*D_map{i}+2.3069;
G0_map{i}=A./F;
end

%$£ Storing parameters and raw data
Other=struct;
Other.Type='Matrix';
Other.ACFx_y=ACFx;
Other.ACFx_y{2}=ACFy;
Other.GSx_y=gs_x;
Other.GSx_y{2}=gs_y;
Other.PHASEx_y=Phase_x;
Other.PHASEx_y{2}=Phase_y;
Other.MODx_y=Mod_x;
Other.MODx_y{2}=Mod_y;
Other.MOD_BKGx_y=Mod_bkg1_x;
Other.MOD_BKGx_y{2}=Mod_bkg1_y;

end
end
function [D,G0,Int,Phase,Modulus]=L_RICS_analysis_Concentration(filename,Mask_LRICS,Mask_Time,step,Range)

if nargin<5
    Range=[1 Inf];
end
   
if ischar(filename)
fprintf(strcat('Loading File...'))
fprintf('\n')
M=LRICS_Line_Import(filename);
else
M=filename;
end

if size(M,1)<size(M,2)
    M=M';
end
fprintf(strcat('Smoothing...'))
fprintf('\n')
M1=LRICS_Line_Smooth_mean0(M,Mask_Time);
M_average=MovAverage_0_line(mean(M),Mask_LRICS);
fprintf(strcat('Analyzing...'))
fprintf('\n')
[Phase,Modulus,Modulus1,~]=LRICS_Line_Concentration(M1(Range(1):min(size(M1,1),Range(2)),:),M_average,Mask_LRICS,step);
L=(Mask_LRICS-1)/2;

Int=mean(M(Range(1):min(size(M1,1),Range(2)),L+1:end-L));
% path=which('LRICS');
% path=path(1:end-7);
path1=which('Param_Cal_6_2_16.mat');
load(path1);
pos=find(mask_vect==L);
D=real(-Param_Cal(pos,3)*log((Phase-Param_Cal(pos,1))/Param_Cal(pos,2)));
F=(Modulus./Modulus1)-1;
A=0.0357*D+2.3069;
G0=A./F;
end
function M1=LRICS_Line_Import(filename)

M=msr2Stack(filename);
[X,Y,Z]=size(M);
M1=permute(M,[2 1 3]);
M1=reshape(M1,[Y,X*Z])';

cnt=1;
while (isempty(nonzeros(M1(cnt,:)))==0)&&(cnt~=X*Z)
    cnt=cnt+1;
end

M1=M1(1:cnt,:);
end
function M=msr2Stack(file_name)

[h_image,~] = omas_bf_open(strcat(file_name));
M = omas_bf_read(h_image,1);
omas_bf_close(h_image);
M=double(M);
if length(size(M))==3
M=permute(M,[2,1,3]);
end
if length(size(M))==4
M=permute(M,[3,2,1,4]);
end
end
function B2=MovAverage_0_line(B,MaskSize)

L=(MaskSize-1)/2;
[X,Y]=size(B);

B2=zeros(1,length(B));
B1=zeros(1,length(B)+2*L);
if X==1
B1(L+1:end-L)=B;
else
B1(L+1:end-L)=B';
end

for i=1:Y
     B2(i)=mean(nonzeros(B1(i:i+2*L)));
end
end 
function M1=LRICS_Line_Smooth_mean0(M,mask)

[X,Y]=size(M);
L=(mask-1)/2;
M1=zeros(X-2*L,Y);

cnt=1;
for i=1+L:X-L
   
    M1(cnt,:)=M(i,:)-mean(M(i-L:i+L,:));
    cnt=cnt+1;
    
end
end        
function [Phase,Modulus,Modulus1,ACF_store]=LRICS_Line_Concentration(M,M_average,mask,step)

[X,Y]=size(M);
if mask==0
    if iseven(Y)
        L=Y/2;
    else
        L=(Y+1)/2;
    end
    Y=L*3;
else
    L=(mask-1)/2;
    Phase=zeros(1,Y-2*L);
    Modulus=zeros(1,Y-2*L);
    Modulus1=zeros(1,Y-2*L);
end
if nargin<4
step=20000;
end
cnt=1;    

for j=1+L:Y-L
    
    if X>step
         n=ceil(X/step);
        for segm=1:n-1
    if mask==0
    [ACF,~]=RICS_ACF2D(M(1+(segm-1)*step:segm*step,:)+M_average(j),L,1,0,0);
    ACF=ACF(L+1,L+1:end);
    else
    [ACF,~]=RICS_ACF2D(M(1+(segm-1)*step:segm*step,j-L:j+L)+M_average(j),L,1,0,0);
    ACF=ACF(L+1,L+1:end);
    end
    if segm==1
        ACF_tmp=ACF;
    else
        ACF_tmp(segm,:)=ACF;
    end
        end
    if mask==0
    [ACF,~]=RICS_ACF2D(M((segm)*step:end,:)+M_average(j),L,1,0,0);
    ACF=ACF(L+1,L+1:end);
    break
    else
    [ACF,~]=RICS_ACF2D(M((segm)*step:end,j-L:j+L)+M_average(j),L,1,0,0);
    ACF=ACF(L+1,L+1:end);
    end
    else
    if mask==0
    [ACF,~]=RICS_ACF2D(M(:,:)+M_average(j),L,1,0,0);
    ACF=ACF(L+1,L+1:end);
    
    gs=fft(ACF);
    gs=conj(gs(2)./abs(gs(1)));
    Phase=atan(imag(gs)./real(gs)); 
    Modulus=sqrt((imag(gs))^2+(real(gs))^2); 
    ACF_store=ACF;
    break
    
    else
        
    [ACF,~]=RICS_ACF2D(M(:,j-L:j+L)+M_average(j),L,1,0,0);
    ACF=ACF(L+1,L+1:end);

    end
    end

    if j==1+L
        ACF_store=ACF;
    else
        ACF_store(cnt,:)=ACF;
    end
    gs=fft(ACF);
    gs=conj(gs(2)./abs(gs(1)));
    
    Phase(cnt)=atan(imag(gs)./real(gs)); 
    Modulus(cnt)=sqrt((imag(gs))^2+(real(gs))^2); 
    gs=fft(ACF+1);
    gs=conj(gs(2)./abs(gs(1)));
    Modulus1(cnt)=sqrt((imag(gs))^2+(real(gs))^2); 
    cnt=cnt+1;
end 
end
function [ACFx_c,ACFy_c,gs_x,gs_y,Phase_x,Phase_y,Mod_x,Mod_y,Mod_bkg1_x,Mod_bkg1_y]=RICS_Local_ACF_MultiImages_Concentration(Matrix,Time_mask,mask,SpatialAverageMask,NumVect)

Matrix=double(Matrix);
M_Average=mean(Matrix,3);
M_Average=MovAverage_0(M_Average,2*mask+1,0);
if Time_mask==0
Matrix=Matrix-repmat(mean(Matrix,3),1,1,size(Matrix,3))+mean(mean(mean(Matrix)));
else
fprintf('Applying Moving Average Time Mask...');
fprintf('\n');
Matrix=RICS_TimeAverage_mean0(Matrix,Time_mask);
end
[X,Y,Z]=size(Matrix);
if NumVect==0
    NumVect=Z;
end
NumVect=NumVect-Time_mask-1;
NumVect(NumVect<=0)=1;

ACFx_c=cell(1);
ACFy_c=cell(1);
A1=zeros(X+2*mask+1,Y,Z);
A1(mask+1:mask+X,1:Y,:)=Matrix;

for i=1:X
    if i==1
        tic
    end
    for j=1:Y-2*mask
        
        A1_tmp=A1(i:i+2*mask,j:j+2*mask,1:NumVect(end))+M_Average(i,j);
        [~,ACF_tmp]=RICS_ACF2D(A1_tmp(:,:,1:NumVect(end)),mask,1,0,1);

        for l=1:length(NumVect)

        % First point is CONSIDERED!!%
        ACFr=mean(ACF_tmp(mask+1,mask+1:end,1:NumVect(l)),3);
        ACFc=mean(ACF_tmp(mask+1:end,mask+1,1:NumVect(l)),3);        
        
        ACFx_c{l}(i,j,:)=ACFr;
        ACFy_c{l}(i,j,:)=ACFc;
        end
        
    end
    if i==1
        T1=toc;
    end
    P=round(T1*(X-i));
    if P>3600
            fprintf(strcat('Processing...',num2str(round((i)/(X)*100)),'/100, Remaining time=',num2str(round(P/3600)),'hour(s)'))
            fprintf('\n')
    else
        if P>60
            fprintf(strcat('Processing...',num2str(round((i)/(X)*100)),'/100, Remaining time=',num2str(round(P/60)),'min'))
            fprintf('\n')
        else
            if P<0
            fprintf(strcat('Processing...',num2str(round((i)/(X)*100)),'/100, Remaining time=0s'))
            fprintf('\n')
            else
            fprintf(strcat('Processing...',num2str(round((i)/(X)*100)),'/100, Remaining time=',num2str(P),'s'))
            fprintf('\n')
            end
        end
    end
end

gs_x=cell(1);
gs_y=cell(1);
gs1_x=cell(1);
gs1_y=cell(1);

for l=1:length(NumVect)
gs_x{l}=fft(ACFx_c{l},[],3);
gs_x{l}=conj(gs_x{l}(:,:,2)./abs(gs_x{l}(:,:,1)));
gs_x{l}(isnan(gs_x{l}))=0;
gs_y{l}=fft(ACFy_c{l},[],3);
gs_y{l}=conj(gs_y{l}(:,:,2)./abs(gs_y{l}(:,:,1)));
gs_y{l}(isnan(gs_y{l}))=0;

gs1_x{l}=fft(ACFx_c{l}+1,[],3);
gs1_x{l}=conj(gs1_x{l}(:,:,2)./abs(gs1_x{l}(:,:,1)));
gs1_x{l}(isnan(gs1_x{l}))=0;
gs1_y{l}=fft(ACFy_c{l}+1,[],3);
gs1_y{l}=conj(gs1_y{l}(:,:,2)./abs(gs1_y{l}(:,:,1)));
gs1_y{l}(isnan(gs1_y{l}))=0;
end

Phase_x{l}=cell(1);
Phase_y{l}=cell(1);
Mod_x{l}=cell(1);
Mod_y{l}=cell(1);
Mod_bkg1_x{l}=cell(1);
Mod_bkg1_y{l}=cell(1);

for l=1:length(NumVect)
Phase_x{l}=atan(imag(gs_x{l})./real(gs_x{l}));
Phase_y{l}=atan(imag(gs_y{l})./real(gs_y{l}));
Mod_x{l}=(sqrt(imag(gs_x{l}).^2+real(gs_x{l}).^2));
Mod_y{l}=(sqrt(imag(gs_y{l}).^2+real(gs_y{l}).^2));
Mod_bkg1_x{l}=(sqrt(imag(gs1_x{l}).^2+real(gs1_x{l}).^2));
Mod_bkg1_y{l}=(sqrt(imag(gs1_y{l}).^2+real(gs1_y{l}).^2));
Phase_x{l}(isnan(Phase_x{l}))=0;
Phase_y{l}(isnan(Phase_y{l}))=0;
Mod_x{l}(isnan(Mod_x{l}))=0;
Mod_y{l}(isnan(Mod_y{l}))=0;    
Mod_bkg1_x{l}(isnan(Mod_bkg1_x{l}))=0;
Mod_bkg1_y{l}(isnan(Mod_bkg1_y{l}))=0;    

Phase_x{l}(real(gs_x{l})<0)=Phase_x{l}(real(gs_x{l})<0)+pi;
Phase_y{l}(real(gs_y{l})<0)=Phase_y{l}(real(gs_y{l})<0)+pi;
end

if SpatialAverageMask~=0
fprintf('Applying Smoothing Mask...');
fprintf('\n');
Phase_y{l}(:,:,2)=MovGaussAverage_sigma(Phase_y{l},SpatialAverageMask*6+1,SpatialAverageMask);
Phase_x{l}(:,:,2)=MovGaussAverage_sigma(Phase_x{l},SpatialAverageMask*6+1,SpatialAverageMask);
Mod_y{l}(:,:,2)=MovGaussAverage_sigma(Mod_y{l},SpatialAverageMask*6+1,SpatialAverageMask);
Mod_x{l}(:,:,2)=MovGaussAverage_sigma(Mod_x{l},SpatialAverageMask*6+1,SpatialAverageMask);
gs_x{l}(:,:,2)=(Mod_x{l}(:,:,2).*cos(Phase_x{l}(:,:,2)))+1i*(Mod_x{l}(:,:,2).*sin(Phase_x{l}(:,:,2)));
gs_y{l}(:,:,2)=(Mod_y{l}(:,:,2).*cos(Phase_y{l}(:,:,2)))+1i*(Mod_y{l}(:,:,2).*sin(Phase_y{l}(:,:,2)));
end
end
function K=MovGaussAverage_sigma(B,MaskSize,sigma)

[X,Y,Z]=size(B);

L=(MaskSize-1)/2;
B1=zeros(X+2*L,Y+2*L,Z);
B1(L+1:end-L,L+1:end-L,:)=B;
for i=1:Z
   tmp=B1(:,:,i);
   tmp(tmp==0)=mean2_0(tmp);
   B1(:,:,i)=tmp;
end

K=zeros(X+2*L,Y+2*L,Z);

for i=1:Z
g=Gaussian2D(MaskSize,sigma);
g=g/sum(sum(g));
K(:,:,i)=conv2(B1(:,:,i),g,'same');
end

K=K(L+1:end-L,L+1:end-L,:);
end
function g2=Gaussian2D(DimMatrix,sigma)

g2=zeros(DimMatrix);
x=(1:1:DimMatrix);
g=Gaussian(x,sigma,(DimMatrix+1)/2,1);

for i=1:1:DimMatrix
    g2(i,1:end)=g*g(i);
end
end
function g=Gaussian(x,sigma,center,height)

g=height*(exp(-(x-center).*(x-center)/(2*sigma.*sigma)));

end
function B=RICS_TimeAverage_mean0(A,Frames)

[X,Y,Z]=size(A);
B=zeros(X,Y,Z-Frames);
cnt=1;

for i=1+(Frames-1)/2:Z-(Frames-1)/2
    
    M=mean(A(:,:,i-(Frames-1)/2:i+(Frames-1)/2),3);
    B(:,:,cnt)=A(:,:,i)-M;
    cnt=cnt+1;

end
end
function K=MovAverage_0(B,MaskSize,Shrink)

K=zeros(1);
[r,c,t]=size(B);
for k=1:t
A=zeros(r+MaskSize-1,c+MaskSize-1);
A((MaskSize-1)/2+1:end-(MaskSize-1)/2,(MaskSize-1)/2+1:end-(MaskSize-1)/2)=B(:,:,k);
for i=1:r
    for j=1:c
        
        H=A(i:i+MaskSize-1,j:j+MaskSize-1);
        M=mean2_0(H);
        K(i,j,k)=M;
        
    end
end

end

if Shrink==1
    K=K((MaskSize-1)/2+1:end-(MaskSize-1)/2,(MaskSize-1)/2+1:end-(MaskSize-1)/2,:);
else
if Shrink==2
    K1=zeros(size(K));
    K1((MaskSize-1)/2+1:end-(MaskSize-1)/2,(MaskSize-1)/2+1:end-(MaskSize-1)/2,:)=K((MaskSize-1)/2+1:end-(MaskSize-1)/2,(MaskSize-1)/2+1:end-(MaskSize-1)/2,:);
    K=K1;
else
end
    
end
end
function A=mean2_0(B)

[~,~,t]=size(B);

A=zeros(1,t);
for i=1:t

A(i)=mean(nonzeros(B(:,:,i)));

end

end
function [Sm,S]=RICS_ACF2D(Sequence,HalfSize,G0corr,Display,Check)

if nargin<4
    Display=0;
    Check=0;
end
if nargin<5
    Check=0;
end

S=ACF2D(Sequence,HalfSize,0);

r=HalfSize+1;
c=HalfSize+1;

if G0corr==1
    S(r,c,:)=S(r,c+1,:);
end

if HalfSize==0
if Check==1
    Sm=zeros(size(S,1),size(S,2));
    count=1;
    for i=1:size(Sequence,3)
        if S(r,c,i)==absmax(S(:,:,i))
            Sm(:,:,count)=S(:,:,i);
            count=count+1;
        end
    end
else
Sm=S;
end
else
    if Check==1
    Sm=zeros(HalfSize*2+1,HalfSize*2+1);
    count=1;
    for i=1:size(Sequence,3)
        if S(r,c,i)==absmax(S(:,:,i))
            Sm(:,:,count)=S(:,:,i);
            count=count+1;
        end
    end
    else
    Sm=S;
    end
end
Sm=mean(Sm,3);
            


if Display==1
    
    Figure_Format_RICS(Sm)
    
%     Smx=Sm(HalfSize+1,HalfSize+1:end);
%     Smy=Sm(HalfSize+1:end,HalfSize+1);
%     Smx=fft(Smx);
%     Smx=Smx(2)/Smx(1);
%     Smy=fft(Smy);
%     Smy=Smy(2)/Smy(1);
    
%    H=HistShow_1(conj([Smx Smy]),100,100,[0,1],[0,1],1);
end
end
function ACF=ACF2D(A,CorrelationRadius,Display)

[X,Y,Z]=size(A);

F=fft2(A);
ACF= F.*conj(F);
ACF=real(fftshift(fftshift(ifft2(ACF),1),2));
G=((sum(sum(A,1),2)).^2/X/Y);

for i=1:Z
    ACF(:,:,i)=ACF(:,:,i)/G(i)-1;
end

if iseven(X)
    r=X/2+1;
else
    r=(X+1)/2;
end
if iseven(Y)
    c=Y/2+1;
else
    c=(Y+1)/2;
end
if CorrelationRadius~=0
    ACF=ACF(r-CorrelationRadius:r+CorrelationRadius,c-CorrelationRadius:c+CorrelationRadius,:);
end

if nargin>3
if Display==1
x=(-CorrelationRadius:1:CorrelationRadius);
figure
surf(x,x,ACF)
zlim([-1,absmax(ACF)])
Figure_Format_Graph
end
end
end
function bool=iseven(x)

if mod(x,2) == 0
bool=1;
else
bool=0;
end
end

function a=absmax(A)

SIZE=ndims(A)-sum(size(A)==1);

switch SIZE
    case 0
        a=A;
    case 1
        a=max(A);
    case 2
        a=max(max(A));
    case 3
        a=max(max(max(A)));
    case 4
        a=max(max(max(max(A))));
end
end
function a=absmin(A)

[r,s,t]=size(A);
SIZE=3;

if t==1
    SIZE=2;
    if r==1||s==1
        SIZE=1;
    end
end

switch SIZE
    case 1
        a=min(A);
    case 2
        a=min(min(A));
    case 3
        a=min(min(min(A)));
end
end




















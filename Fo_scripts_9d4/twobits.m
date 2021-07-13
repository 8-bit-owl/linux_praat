files = ['V:\Offsite Data Collection\Lakeshore Data\exp2 preped\to run 44100\est'];
%code = ['V:\1 General Analysis And Techniques\Matlab Universal Scripts\Fo_scripts'];
%code = ['V:/1 General Analysis And Techniques/Matlab Universal Scripts/Fo_scripts'];
%code = ['/mnt/ufs18/home-070/ganjaisa/Acoustics/new_praat/Git_Praat']
code = ['./']
addpath(code)
cd(files) 
filename = FilenamesByExt('wav');
start = 1;
cnt=length(filename);
for j_iter = start:cnt  
    fname=filename{j_iter};
    [y,Fs] = audioread(fname);
    LTAS=ltas(y,Fs,min([2^floor(log2(length(y))) 2^13]),'han',0,0,1,-100);
        oct3 = LTAS.oct3lev;
        f = (1000).*((2^(1/3)).^[-13:13]);
        octs = (1:27)/3;
    p1 = polyfit(f,oct3',1);
    LTAS.tiltf = p1(1);
    p2 = polyfit(octs,oct3',1);
    LTAS.tilto = p2(1);
    out{1,j_iter} = fname;
    out{2,j_iter} = LTAS.tiltf;
    out{3,j_iter} = LTAS.tilto;
end
    
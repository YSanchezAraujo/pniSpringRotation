using Plots;
using Distributions, Wavelets, CMF, DataFrames, MAT, DSP, PyPlot;
using PyCall;

scipySig = pyimport("scipy.signal");

utilsPath = "/home/yoel/Desktop/pniSpringRotation/scripts/utils.jl"
include(utilsPath);

dataPath = "/home/yoel/Desktop/pniSpringRotation/data/data_20150120F1";
# information on the electrodes, to be used for plotting eventually
eleInfo = joinpath(dataPath, "ElectrodesFr.mat");

# full path to data
dataFiles = [joinpath(dataPath, x) for x in readdir(dataPath) if occursin("ECoG_ch", x)];

# string parser function for the load_all_data function
parseFn(x::String) = parse(Int64, split(split(split(x, "/")[end], ".")[1], "ch")[end])

# load in data
data = load_all_data(dataFiles, parseFn);

FS = 1000;
w = 10.;
freqs = collect(1:1:(FS/20));
widths = w*FS ./ (2*freqs*pi);
cwtm = scipySig.cwt(data[:, 14], scipySig.morlet2, widths, w=w);
scwtm = abs.(cwtm);

# downsample
ratio = 1//10;
rscwtm = downsample_spectro(scwtm, size(data, 1), ratio);

# one time plot of the entire thing
show_spectro(scwtm, (15, 8), "test")

# for making a gif
ranges = 1:1000:size(rscwtm, 2);
nr = length(ranges);

anim = @animate for r in 2:nr
    samp_idxs = ranges[r-1]:ranges[r]
    heatmap(rscwtm[:, samp_idxs], xlabel="samples (FS=1000Hz)", ylabel="Freq (Hz)")
end

gif(anim, "test2.gif", fps=2)


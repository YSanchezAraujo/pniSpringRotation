using Plots;
using Distributions, Wavelets, CMF, DataFrames, MAT, DSP, PyPlot;
utilsPath = "/home/yoel/Desktop/pniSpringRotation/scripts/utils.jl"
include(utilsPath);

DATA_SAMPLING_RATE = 1000; # Hz
dataPath = "/home/yoel/Desktop/pniSpringRotation/data/data_20150120F1";
# information on the electrodes, to be used for plotting eventually
eleInfo = joinpath(dataPath, "ElectrodesFr.mat");

# full path to data
dataFiles = [joinpath(dataPath, x) for x in readdir(dataPath) if occursin("ECoG_ch", x)];

# string parser function for the load_all_data function
parseFn(x::String) = parse(Int64, split(split(split(x, "/")[end], ".")[1], "ch")[end])

# load in data
data = load_all_data(dataFiles, parseFn);
NSAMPLES = size(data, 1);
NSECONDS = NSAMPLES / DATA_SAMPLING_RATE;

# create a morlet wavelet
# the number inside WT.Morlet() is related to the trade-off between
# resolutions, the smaller it is, the less we pick up from the lower frequencies
# at about 3, you get a small amount from around the 10 Hz range
# the parameter s, seems to control how many frequencies are computed
# at s=20 we get up to 71hz, lower a smaller s it's less, at s = 30 it's 108hz
wav = wavelet(WT.Morlet(6); s=30);
# example of computing a continuous time morlet wavelet transform on an array
cwtMor = cwt(data[:, 1], wav);
# extracting the real part of the transform
rcwtMor = real(cwtMor);

# for downsampling
ratio = 1//10;
rsTransform = downsample_spectro(rcwtMor, size(data, 1), ratio);

for ch_num in 1:size(data, 2)
    cwtMor = cwt(data[:, ch_num], wav)
    rcwtMor = real(cwtMor)
    rscwt = downsample_spectro(rcwtMor, size(data, 1), ratio)
    save_name = string("fullspectro_", ch_num, ".png")
    show_spectro(abs.(rcwtMor'), (15, 8), split(save_name, ".")[1])
    plt.savefig(save_name, dpi=100, bbox_inches="tight")
    plt.close()
    save_name = string("dsspectro_", ch_num, ".png")
    show_spectro(abs.(rscwt'), (15, 8), split(save_name, ".")[1])
    plt.savefig(save_name, dpi=100, bbox_inches="tight")
    plt.close()
end

# remember to transpose at plot time so that time is on the x-axis
# creating a gif of the spectrogram

ranges = 1:1000:size(rsTransform, 1);
nr = length(ranges);

anim = @animate for r in 2:nr
    samp_idxs = ranges[r-1]:ranges[r]
    heatmap(abs.(rsTransform'[:, samp_idxs]), xlabel="time (s)", ylabel="Freq (Hz)")
end

gif(anim, "test2.gif", fps=2)

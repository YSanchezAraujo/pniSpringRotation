using JLD, Distributions, Wavelets, CMF, DataFrames, MAT;
utilsPath = "/home/yoel/Desktop/pniSpringRotation/scripts/utils.jl"
include(utilsPath);

dataPath = "/home/yoel/Desktop/pniSpringRotation/data/data_20150120F1";
eleInfo = joinpath(dataPath, "ElectrodesFr.mat");

# full path to data
dataFiles = [joinpath(dataPath, x) for x in readdir(dataPath) if occursin("ECoG_ch", x)];

# string parser function for the load_all_data function
parseFn(x::String) = parse(Int64, split(split(split(x, "/")[end], ".")[1], "ch")[end])

# load in data
data = load_all_data(dataFiles, parseFn);

# create a morlet wavelet struct
morlet = Wavelets.WT.Morlet(3);

# example of computing a continuous time morlet wavelet transform on an array
cwtMor = cwt(data[:, 1], morlet);

# this computes it for 32-34 different frequencies it seems, if I'm understanding this correctly

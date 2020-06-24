using Plots;
using Distributions, CMF, DataFrames, MAT, PyPlot;
using PyCall;
using DSP;

scipySig = pyimport("scipy.signal");

utilsPath = "/home/yoel/Desktop/pniSpringRotation/scripts/utils.jl"
include(utilsPath);

dataPath = "/home/yoel/Desktop/pniSpringRotation/data/data_20150120F1";
# information on the electrodes, to be used for plotting eventually
eleInfo = joinpath(dataPath, "ElectrodesFr.mat");
eventsPath = joinpath(dataPath, "repEvent.mat");
events = matread(eventsPath);

# full path to data
dataFiles = [joinpath(dataPath, x) for x in readdir(dataPath) if occursin("ECoG_ch", x)];

# string parser function for the load_all_data function
parseFn(x::String) = parse(Int64, split(split(split(x, "/")[end], ".")[1], "ch")[end])

# load in data
data = load_all_data(dataFiles, parseFn);

FS = 1000;
w = 10.;
nA = 30;

avgs = zeros(201, 31);
colidx = setdiff(1:32, 18);
for (k, pidx) in enumerate(colidx)
    avgs[:, k] = avg_event(5.0, data[:, pidx], events["repEvent"])
end

fig, ax = plt.subplots(ncols=1, nrows=3, figsize=(15, 12));
for (n, m) in enumerate([3.0, 5.0, 11.0])
    avgs = zeros(201, 31);
    for (k, pidx) in enumerate(colidx)
	avgs[:, k] = avg_event(m, data[:, pidx], events["repEvent"])
    end
    ax[n].plot(avgs, alpha=0.5, lw=0.9)
    ax[n].plot(mean(avgs, dims=2), color="black", lw=2.0)
end

#concated_scwtm = zeros(33*size(data, 2), size(data, 1));
#
#c = 1;
#for k in 1:size(data, 2)
#    scwtm = get_spectro(data[:, k], w, FS, nA);
#    concated_scwtm[c:(c+32), :] = scwtm
#    global c += 33
#end

scwtm = get_spectro(data[:, 1], w, FS, nA);
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
    Plots.plot(
	 heatmap(rscwtm[:, samp_idxs], xlabel="samples (FS=1000Hz)", ylabel="Freq (Hz)"),
	 heatmap(rsv[:, samp_idxs], xlabel="samples (FS=1000Hz)", ylabel="Freq (Hz)")
    )
end

gif(anim, "recontest2.gif", fps=2)

meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))

# fitting cNMF

res = fit_cnmf(rscwtm; L=100, K=8, alg=HALSUpdate, l1_H=0.4, l2_H=0.2, l1_W=0.4, l2_W=0.5);


function make3drecon(W, H, k)
    K, F, L = size(W);
    T = size(H, 2);
    lrecons = zeros(F, L, T);
    for l in 1:L
	lrecons[:, l, :] = W[k, :, l]*H[k, :]'
    end
    return lrecons
end

test = make3drecon(res.W, res.H, 1);
# for surface plotting
np = pyimport("numpy");
z, y, x = np.nonzero(np.array(test[:, :, 1:1000]));

fig = plt.figure();
x, y = meshgrid(1:33, 1:100);
i=300;
for k in [311, 312, 313]
    ax = fig.add_subplot(k, projection="3d");
    ax.plot_surface(reshape(x, (33, 100)), reshape(y, (33, 100)), test[:, :, i], cmap="viridis");
    ax.set_zlabel("power");
    ax.set_xlabel("freq (Hz)");
    ax.set_ylabel("L length");
end
plt.savefig("K3RECON.png", dpi=100, bbox_inches="tight");
ax.plot_trisurf(x, z, y, cmap="viridis", edgecolor="none");
ax.set_zlabel("power");
ax.set_xlabel("freq (Hz)");
ax.set_ylabel("L length");
plt.savefig("K3_recon.png", dpi=100, bbox_inches="tight");
rsv = downsample_spectro(v, size(data, 1) ratio);

#fig, ax = plt.subplots(nrows=1, ncols=5, figsize=(15, 9))
#for k in 1:5
#    im = ax[k].imshow(results.W[k, :, :])
#    ax[k].set_title(string("Component ", k))
#    fig.colorbar(im, ax=ax[k], fraction=0.03)
#end
#
#
#c1proj = results.W[1, :, 1] * results.H[1, :]'
#rsc1proj = downsample_spectro(c1proj, size(data, 1), ratio);
#
#anim = @animate for r in 2:nr
#    samp_idxs = ranges[r-1]:ranges[r]
#    heatmap(rsc1proj[:, samp_idxs], xlabel="samples (FS=1000Hz)", ylabel="Freq (Hz)")
#end
#gif(anim, "recontestc1.gif", fps=2)


function load_all_data(paths::Array{String, 1}, f::Function)
    nFiles = length(paths)
    # load in first file to get number of samples
    data = vcat(collect(values(matread(paths[1])))...)'
    # get the ecog array number
    arrNum = f(paths[1])
    # for sorting columns later
    arrOrder = [arrNum]
    
    for k in 2:nFiles
	# get substring to order the columns by ecog array number
	arrNum = f(paths[k])
	# save order
	push!(arrOrder, arrNum)
	# load in data
	vars = vcat(collect(values(matread(paths[k])))...)'
        # concatenate the data
	data = hcat(data, vars)
    end

    return data[:, sortperm(arrOrder)]
end

# this works sorta ok
function viz_ecog(info_path::String; interval=5)
    fr = matread(info_path)
    outline = Int64.(fr["LINE"]);
    outline = (outline .- 255) ./ 255
    X, Y = Int64.(fr["X"]), Int64.(fr["Y"])
    
    for (idx, pair) in enumerate(zip(Y, X))
	p1int = collect((pair[1] - interval):(pair[1] + interval))
	p2int = collect((pair[2] - interval):(pair[2] + interval))
	for pp in zip(p1int, p2int)
	    outline[pp[1], pp[2], :] .= 1.0
	end
    end

    return outline
end

function show_spectro(data, figsize, title)
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.matshow(data, origin="lower")
    ax.axis("tight")
    ax.xaxis.set_ticks_position("bottom")
    ax.set_ylabel("Frq (Hz)")
    ax.set_xlabel("Samples (Fs=1000Hz)")
    ax.set_title(title)
    fig.colorbar(im)
end

function downsample_spectro(wav_res, true_N, ratio)
    nfreqs = size(wav_res, 1)
    nSamples = Int64(ceil(Float64(true_N * ratio)))
    rsTransform = zeros(nfreqs, nSamples);
    for k in 1:nfreqs
        rsTransform[k, :] = resample(wav_res[k, :], ratio)
    end
    return rsTransform
end


#
#function save_giff(data, saveName, fps)
#    @userplot rwPlot
#    @recipe function g(gs::rwPlot)
#        i = gs.args
#        plot(i)
#    end
#
#    anim = @animate for i in data
#        plot(i)
#    end
#
#    gif(anim, saveName, fps=fps)
#end
#

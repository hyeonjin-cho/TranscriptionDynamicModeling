using DelimitedFiles
using Distributed
using Pkg
#Pkg.add(name="StochasticGene", version="1.1.7")

@everywhere using StochasticGene
Pkg.status("StochasticGene") # to confirm v1.1.7

gene = ARGS[1]
cell = "HBEC"
root = "."
decayrate = StochasticGene.get_decay(gene,cell,root)
transitions = ([1,2], [2,1], [2,3], [3,1]) # KP model
G = 3
R = 3
S = R
insertstep = 1
root = "."
splicetype = ""
datatype = "rnadwelltime"
onstates = [Int[], Int[], [2,3], [2,3]]
dttype=["ON", "OFF", "ONG", "OFFG"]
nhist = 600
nalleles = StochasticGene.alleles(gene, cell, root)
propcv = 0.3 # adjust this value to get the optimal acceptance rate
priorcv = 10.0
cv = 10.0
nchains = 10 # adjust this value to get the optimal rhat value (<1.1)
samplesteps = 20000
warmupsteps = 10000
annealsteps = 0
maxtime = 20000.
temp = 1.
tempanneal = 1.
tempfish = 1
type = ""
optimized = 0
bs = 0
writesamples = true
label = "gt"
rtype = "ml"
method = 1
options = StochasticGene.MHOptions(samplesteps,warmupsteps,annealsteps,maxtime,temp,tempanneal)
resultfolder = string(gene,"_results")
resultfolder = StochasticGene.folder_path(resultfolder, root, "results", make=true)
counts = 1
total=1000000
tol=1e-6
hierarchical = tuple()
noisepriors = []
noiseparams = occursin("trace", lowercase(datatype)) ? length(noisepriors) : zero(Int)
#priormean = StochasticGene.prior_ratemean(transitions, R, S, insertstep, decayrate, noisepriors)
priormean = [5, 5, 5, 5, 0.1, 0.1, 0.1, 0.1, 0.05, 0.05, 0.05, decayrate]
fittedparam = [1,2,3,4,5,8,9]
fixedeffects = ([5,6,7], [9,10,11])
probfn = StochasticGene.prob_GaussianMixture
weightind = 5

#=================Load data and model=================#
infile = string("data/HBEC_offTime/"*gene*"_ON_OFF_hist_PDF.csv")
LC = readdlm(infile,',')
x = StochasticGene.truncate_histogram(LC[:,2],.999,1000)
LC = LC[1:length(x),:]
infile = string("data/HBEC_SC/PiC_TBP_survivalparams.csv")
SC = readdlm(infile,',')[1:493,1:2]
SC = convert(Matrix{Float64}, SC)
x = StochasticGene.truncate_histogram(SC[:,2],0.999,1000)
SC = SC[1:length(x),:]
fishfile = string("data/HBEC_smFISH/"*gene*"_steady_state_mRNA.csv")
x = readdlm(fishfile, ',')[3:end,3]
x = x[x .!= ""]
x = StochasticGene.truncate_histogram(x,.99,1000)
x /= sum(x)
x *= counts
histFISH = x

#=================Adjust LC and SC data=================#
tempLC = zeros(Float64, 2000, 3)
LC = vcat(LC, tempLC)
bins_LC = collect(5/3:5/3:size(LC)[1]*5/3)
tempSC = zeros(Float64, 2000, 2)
SC = vcat(SC, tempSC)
bins_SC = collect(7/300:1/300:size(SC)[1]*1/300+20/900)
bins = [bins_LC, bins_LC, bins_SC, bins_SC]
h = [LC[:,2], LC[:,3], SC[:,2], SC[:,2]]

name = string("_",label,"_",gene,"_",G,R,S,insertstep,"_",nalleles,".txt")
r = StochasticGene.readrates(joinpath(resultfolder, "rates" * name), 1)

if isempty(r)
    r = priormean
end
println(gene)
println(r)

data = StochasticGene.RNADwellTimeData(label, gene, length(histFISH), histFISH, bins, h, dttype)
model = StochasticGene.load_model(data, r, priormean, fittedparam, fixedeffects, transitions, G, R, S, insertstep, nalleles, priorcv, onstates, decayrate, propcv, splicetype, probfn, noiseparams, weightind, hierarchical)
StochasticGene.print_ll(data,model)

fittt, stats, measures = StochasticGene.run_mh(data, model, options, nchains)
StochasticGene.finalize(data,model,fittt,stats,measures,options.temp,resultfolder,optimized,bs,writesamples)

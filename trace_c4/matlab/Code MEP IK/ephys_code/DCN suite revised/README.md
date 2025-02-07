# CSpk Suite

## Detection of switch sessions, keeping into account both pure behavior and ephys sessions

look at:
    "/home/narain/BayesLab Dropbox/Julius Koppen/TraceExperiments/BehavioralAnalysis/BlinkAnalysis_Behavior_AcrossMice_IntegrationwithEphys_III.m"

## pre-Curation

- extract simple spikes :  experiment output data (raw) -> Simple Spike idcs

- extract complex spike candidates (k-means clustering + base rate filtering) : raw -> Complex Spike Candidates

- cross-correlation : CSpk candidates -> cross-corr results + cross-cor plots (for interactive scoring)

- contamination analysis : CSpk candidates -> contamination values (inspects Simple Spike pause)

- CSpk PSTH baseline and SD filtering : CSpk candidates -> Bool

- CSpk occurrence percentage (across trials) filtering : CSpk candidates -> Bool

- final verdict: {cross-corr, PSTH filter, contamination score, manual inspection} -> Bool

## Curation


## Figures

- histograms of neurons pooled together for all mice within condition
	- all Delta Expert vs Uniform Expert
	- all DEUN Delta Expert vs DEUN Uniform Novice
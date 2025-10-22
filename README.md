**Model Program Versions**:

LHC-Walk-v01.R: No clustering, only heatmap

LHC-Walk-v02.R: Migration is an event not piecewise. Saves individual graphs for each simulation

LHC-Walk-v03.R: Added R0 calculation (NGM)

LHC-Walk-v04.R: First attempt at Heirarchical Euclidean Clustering

LHC-Walk-v05.R: Separated Albatross and Predator migration. Second attempt at Heirarchical Euclidean Clustering

LHC-Walk-v06.R: Got the validation functions up and running (most importantly validate_Births) and got them working correctly

LHC-Walk-v07.R: Implemented here() library for proper path relativization. Fixed timings to be based on literature. Fixed births to be realistic. Changed seasonality
                curves to match Dr. Gamble's description. Feature extraction. Final attempt at Euclidean Clustering.

LHC-Walk-v08.R: Histogram-based K-means attempt 1 (see Advanced Clustering Scripts for furhter work done on these data sets)

**Advanced Clustering Scripts:**

_Note: All Advanced Clustering Scripts used data from LHC-Walk-v08.R_

ClusteringExperiment.R - K-means feature-based clustering (just peak curves for infection)

Clustering_Experiment-2.R - K-means on species-relativized curves for Infected and Dead compartments of all three species

Clustering_Experiment-3.R - "Clustering" by divergence of features stdevnum away from mean for all samples

Clustering_Experiment-4.R - K-means clustering by features (with some excluded because of their relative similarity)

Clustering_Experiment-5.R - K-means using KL distance rather than WCSS as optimization criterion.

Clustering_Experiment-6.R - Dirichlet Process Clustering, cannot be run on local machine, needs HPC validation beofre use.

VariableExploration-v01/2.R - Used to determine what variables are too similar for Clustering_Experiment-4.R

**Scenario Analysis Scripts:**

_Note: All Scenario Analysis Scripts used data from LHC-Walk-v08.R_


Scenario 1: Overwintering

Scenario 2: Multiple peaks

# Package PeptideClassifier

## Description

The main purpose of this package is to facilitate the exploration of sets of peptides using a flexible shiny app. The oligopeptides can be grouped using various clustering techniques.

**Steps:**
- Start Shiny application: runApp();
- Import oligo-peptides from an external fasta file;
- Browse and cluster the peptides;
- Visualise the clusters;

The peptides are split into various n-Grams. Under the **DTM** tab are the controls to construct a DocumentTermMatrix, which will be used during the subsequent steps.

Topic modelling using LDA is available under the **Topics** tab. Peptides binding to the same site on the same receptor would correspond to a specific "topic". Correlated Topic Modeling may be useful to detect more complicated binding sites or homologous receptors with slight variation of the binding site.

Hierarchical clustering can be performed under the **Clustering** tab. Previously saved clusters/trees can be loaded into the app, thus avoiding expensive computation times.

The **Diagnostics** and **Correlation** tabs enable the comparison of the various hierarchical clustering algorithms. Warning: Correlation between trees takes hours to compute!


### Rational

The main purpose of this tool is to group peptides together, not to predict a certain bioactivity of the peptides.

The main rational is as follows: peptides (with a specific bioactivity) will bind to certain receptors on the surface of the cells. There is a finite number of receptors (e.g. receptors active in inflammation, see example below). Each receptor has probably a limited number of sites where inhibitors can block efficiently that receptor. The total number of binding sites on these receptors is therefore a finite number.

Peptides binding to the same site probably share some structural/chemical features which facilitate that binding. The purpose of this tool is to group peptides with specific bio-activities into a number of clusters, based on the assumption that peptides binding around a common site share common structural/chemical features. Each cluster would correspond to such a binding site - although the number of such binding sites is not known a-priori.

Note: peptides binding to a receptor may be oriented in various ways. The user can select both **directed** and **undirected** n-Grams: the undirected n-Grams enable the simulation of binding of a peptide rotated by 180 degrees. The assumption is that the peptide backbone plays only a minor role in the binding to a receptor; the side chains of the various amino-acids play probably the biggest part and the order of these amino acids is reversed if the peptide is rotated by 180 degrees.


### Chemical Descriptors

A very limited set of chemical descriptors is currently implemented. The plan is to extend the number of these descriptors - as time permits (e.g. with hydrophobicity indexes).

Note: the chemical descriptors are calculated on n-Grams of specific sizes. These descriptors are per-se undirected.



### Examples

A dataset with anti-inflammatory peptides (AIP) is available in the /inst/examples folder. A number of precomputed trees is also available in this folder.

Note:
- The dataset was copied from the BertAIP project on GitHub. The description for that project states that BertAIP is a bidirectional encoder representation from transformers-based tool for the prediction of anti-inflammatory peptides (see also Ref below).
> https://github.com/ying-jc/BertAIP
- The training and testing datasets have been joined together (see file AIP.all_dataset.fasta).
- The original training dataset is also available on its own (file AIP.train_dataset.fasta); but note that this application is **not** used for such predictions.
- Although epitopes bind directionally in the groove of the MHC I & MHC II molecules, the corresponding TCRs may each be oriented differently.
- MHC-Epitope-TCR Complexes: may be a very specialised example, which is less suitable for the more general purpose of this application. But the datasets were readily available and fit my interests in immunology.

The construction of the dataset is described in:
- Xu T, Wang Q, Yang Z, Ying J. A BERT-based approach for identifying anti-inflammatory peptides using sequence information. Heliyon. 2024 Jun 13;10(12):e32951. PMID: 38988537.
> https://doi.org/10.1016/j.heliyon.2024.e32951.


## Authors

1. PeptideClassifier: Author, Maintainer: L. Mada


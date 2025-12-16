# run redirecting to stdout to generate GO file
# --> python GO
# _transcription_factors.py > "output_GO_list_TFs.txt"

import pandas as pd
from local_utils import get_gene_ontology
import sys

data = pd.read_csv("whole_Muscle-regulons.csv", sep='\t', header=0) # Dataset from SCENIC (https://scenic.aertslab.org)

#Finding all TFs that regulate each of these targets
tfs_Luc7l3 = list(data.TF[data.gene == "Luc7l3"])
tfs_Cd59a = list(data.TF[data.gene == "Cd59a"])
tfs_Fabp3 = list(data.TF[data.gene == "Fabp3"])
tfs_Pcnp = list(data.TF[data.gene == "Pcnp"])

tfs_all = tfs_Luc7l3 + tfs_Cd59a + tfs_Fabp3 + tfs_Pcnp
tfs = set(tfs_all)

get_gene_ontology(tfs, species="mouse")
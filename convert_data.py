import pandas as pd
import random
import matplotlib.pyplot as plt
import seaborn as sns

def jaccard(l1, l2):
    return len(l1.intersection(l2)) / len(l1.union(l2))

tfs = []
genes = []
with open("c3.tft.v2023.2.Hs.symbols.gmt", "r") as file:
    line = file.readline().strip()
    while(line != ""):
        line = line.split('\t')
        tf = line[0].split("_")[0]
        gene_per_tf = line[2:]
        gene_per_tf = [i.upper() for i in gene_per_tf]
        tfs.extend([tf]*len(gene_per_tf))
        genes.extend(gene_per_tf)
        line = file.readline().strip()

tf_df = pd.DataFrame()
tf_df["TF"] = tfs
tf_df["GENEID"] = genes
tf_df.to_csv("tf_gene_mapping.csv", index = False)

gos = []
genes = []
with open("c5.go.bp.v2023.2.Hs.symbols.gmt", "r") as file:
    line = file.readline().strip()
    while(line != ""):
        line = line.split('\t')
        go = line[0]
        gene_per_go = line[2:]
        gene_per_go = [i.upper() for i in gene_per_go]
        gos.extend([go]*len(gene_per_go))
        genes.extend(gene_per_go)
        line = file.readline().strip()

go_df = pd.DataFrame()
go_df["GO"] = gos
go_df["GENEID"] = genes
go_df.to_csv("go_gene_mapping.csv", index = False)

merged_df = pd.merge(tf_df, go_df, how = "inner", on = "GENEID").reset_index()[["TF", "GENEID", "GO"]]
merged_df.to_csv("merged_df.csv", index = False)

tf_go = merged_df[["TF", "GO"]].copy().drop_duplicates().reset_index()[["TF", "GO"]]
tfs_only = merged_df[["TF"]].copy().drop_duplicates().reset_index()[["TF"]]
tfs_only = tfs_only.merge(tfs_only, how = "cross")
tfs_only = tfs_only[tfs_only["TF_x"] != tfs_only["TF_y"]].reset_index()[["TF_x", "TF_y"]]

tf_xs = []
tf_ys = []
jaccard_distances = []
for i in range(40000):
    rand_iloc = random.randrange(len(tfs_only))
    print(i)
    tfx_gos = set(tf_go[tf_go["TF"] == str(tfs_only.iloc[rand_iloc]["TF_x"])]["GO"].values.tolist())
    tfy_gos = set(tf_go[tf_go["TF"] == str(tfs_only.iloc[rand_iloc]["TF_y"])]["GO"].values.tolist())
    jaccard_distances.append(jaccard(tfx_gos, tfy_gos))
    tf_xs.append(tfs_only.iloc[rand_iloc]["TF_x"])
    tf_ys.append(tfs_only.iloc[rand_iloc]["TF_y"])

tfs_only_subset = pd.DataFrame()
tfs_only_subset["TF_x"] = tf_xs
tfs_only_subset["TF_y"] = tf_ys
tfs_only_subset["JACCARD"] = jaccard_distances
print(tfs_only_subset["JACCARD"].describe())
tfs_only_subset.to_csv("jaccard_permtest.csv", index = False)

s4 = pd.read_csv("./s4.csv").set_index("GO")
colnames = list(s4.columns)
tf_xs = []
tf_ys = []
jaccards = []

for i in range(len(colnames)):
    for j in range(len(colnames)):
        tf_xs.append(colnames[i])
        tf_ys.append(colnames[j])
        jaccards.append(jaccard(set(list(s4[s4[colnames[i]] == 1].index)), set(list(s4[s4[colnames[j]] == 1].index))))

pairwise_df = pd.DataFrame()
pairwise_df["TF_x"] = tf_xs
pairwise_df["TF_y"] = tf_ys
pairwise_df["Jaccard"] = jaccards
pairwise_df.to_csv("./pairwise_jaccard.csv")

pairwise_jaccard = pairwise_df.pivot(index = "TF_x", columns = "TF_y", values = "Jaccard")
sns.clustermap(pairwise_jaccard, linewidth = 0.5, cmap = "viridis")
plt.savefig("./clustermap.pdf")
#plt.show()

#comparison_df = pairwise_df.unstack().reset_index().rename(columns = {"level_0": "TF_x", "TF": "TF_y", 0: "Jaccard"})
#comparison_df = comparison_df[comparison_df["Jaccard"] != 1]
#comparison_df = comparison_df.reset_index()[["TF_x", "TF_y", "Jaccard"]]
#comparison_df["pval"] = [len(jaccard_permtest[jaccard_permtest["JACCARD"] > comparison_df.iloc[i]["Jaccard"]])/5000 for i in range(len(comparison_df))]
#comparison_df.sort_values(by = "pval", inplace = True)
#comparison_df = comparison_df.pivot(index = "TF_x", columns = "TF_y", values = "pval").fillna(1)

pairwise_df["pval"] = [len(tfs_only_subset[tfs_only_subset["JACCARD"] > pairwise_df.iloc[i]["Jaccard"]])/5000 for i in range(len(pairwise_df))]
#pairwise_df["padj"] = pairwise_df["pval"] / ((len(colnames)*len(colnames) - len(colnames))/2)
pvals = []
for i in range(len(pairwise_df)):
    pvals.append(1 if pairwise_df.iloc[i]["TF_x"] == pairwise_df.iloc[i]["TF_y"] else pairwise_df.iloc[i]["pval"]/((len(colnames)*len(colnames) - len(colnames))/2))

pairwise_df["padj"] = pvals
pairwise_df.to_csv("./pairwise_jaccard_post_correction.csv")

pairwise_df_pval = pairwise_df.pivot(index = "TF_x", columns = "TF_y", values = "padj")
sns.heatmap(pairwise_df_pval, cmap = ["red", "blue"], center = 0.05/105, linewidth = 0.5)
plt.savefig("./pvalmap.pdf")
#plt.show()

sns.histplot(tfs_only_subset["JACCARD"])
plt.savefig("./jaccard_permtest_histplot.pdf")
#plt.show()

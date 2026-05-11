import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

### set of the arguments for the script, the first one is the path of the input file,
### the second one is the path of the output directory, where the results will be saved

input_path = sys.argv[1]
output_path = sys.argv[2]
nomefile = os.path.basename(input_path)
nomebase = os.path.splitext(nomefile)[0]

df = pd.read_csv(input_path, sep=r"\s+", header=None, names=["Time", "Potential", "Total_energy", "Temperature", "Pressure"], na_values=["NaN", "nan", "NA", "N/A", "None", "null", "NULL", ""])
missing = df.isna().sum()
### the df.isna() function returns a boolean DataFrame where True indicates missing values, 
### and the sum() function counts the number of missing values in each column. 
### The resulting Series 'missing' contains the count of missing values for each column
### in the DataFrame.
### in this case was used as a mild control 
if missing.any():
    print("Error: invalid values found.")
    print(missing[missing > 0])
    raise SystemExit(1)
print(df.isna().any().any())

headers_full = df.columns.tolist()
headers = headers_full[1:]

risultati = {}

for header in headers:
    array = df[header].to_numpy()

    risultati[header] = {
        "Media": np.mean(array),
        "Massimo": np.max(array),
        "Minimo": np.min(array),
        ### ddof = 1 is used to calculate the sample standard deviation, 
        ### which divides by (n - 1) instead of n, where n is the number of observations.
        ### This method it was chosen because it's a little bit more conservative and it 
        ### gives a better estimate of the population standard deviation when the sample 
        ### size is small, which is the case of our simulations.
        "Deviazione_standard": np.std(array, ddof=1)
    }

tabella_risultati = pd.DataFrame.from_dict(risultati, orient="index").T


output_file = os.path.join(output_path, f"{nomebase}_tabella.tsv")
tabella_risultati.to_csv(output_file, sep="\t", index=True)



fig, ax = plt.subplots(2,2, figsize=(12, 8))
fig.subplots_adjust(wspace=0.5, hspace=0.6, right=0.8)

for asse in ax.flat:
    asse.grid()
    asse.set_xlabel("time (ps)", fontsize=8)
    asse.tick_params(labelsize=8)


ax[0,0].plot(df["Time"], df["Temperature"], label="Temperature", color="blue", alpha=0.5)
ax[0,0].set_ylabel("Temperature (k)", fontsize=8)
ax[0,1].plot(df["Time"], df["Total_energy"], label="Total_energy", color="green", alpha=0.5)
ax[0,1].set_ylabel("Total Energy (kj/mol)", fontsize=8)
ax[1,0].plot(df["Time"], df["Pressure"], label="Pressure", color="red", alpha=0.5)
ax[1,0].set_ylabel("Pressure (bar)", fontsize=8)
ax[1,1].plot(df["Time"], df["Potential"], label="Potential", color="orange", alpha=0.5)
ax[1,1].set_ylabel("Potential (kj/mol)", fontsize=8)
fig.legend(loc="center left", bbox_to_anchor=(0.85, 0.9), fontsize=8, frameon=False)
output_file = os.path.join(output_path, f"{nomebase}_grafici.png")
fig.savefig(output_file)
plt.close(fig)

tabella_risultati_round = tabella_risultati.round(2)

fig, ax = plt.subplots(figsize=(17, 5))
ax.axis("off")

tabella = ax.table(
    cellText=tabella_risultati_round.values,
    rowLabels=tabella_risultati_round.index,
    colLabels=tabella_risultati_round.columns,
    loc="center",
    colColours=["lightblue"] * len(tabella_risultati_round.columns),
    rowColours=["lightgray"] * len(tabella_risultati_round.index)
    
)

tabella.auto_set_font_size(False)
tabella.set_fontsize(10)
tabella.scale(1, 1.5)
output_file = os.path.join(output_path, f"{nomebase}_imgtab.png")
fig.savefig(output_file)
plt.close(fig)

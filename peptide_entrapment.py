#  Licensed to the Apache Software Foundation (ASF) under one
#  or more contributor license agreements.  See the NOTICE file
#  distributed with this work for additional information
#  regarding copyright ownership.  The ASF licenses this file
#  to you under the Apache License, Version 2.0 (the
#  "License"); you may not use this file except in compliance
#  with the License.  You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing,
#  software distributed under the License is distributed on an
#  "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
#  KIND, either express or implied.  See the License for the
#  specific language governing permissions and limitations
#  under the License.

import glob
import os
import pandas as pd


def main(args):
    entrapment_prefix = args[0]
    wd = args[1]
    tool_name = args[2]
    r = float(args[3])

    target_peptide_count, entrapment_peptide_count = calculate(wd, entrapment_prefix, tool_name)

    print(tool_name + ":")
    print("Global-wise peptides:")
    print(f"Target: {target_peptide_count}")
    print(f"Entrapment: {entrapment_peptide_count}")
    print(f"combined method: {entrapment_peptide_count * (1 + 1 / r) * 100 / (target_peptide_count + entrapment_peptide_count):.2f}%")
    print(f"lower bound: {entrapment_peptide_count * 100 / (target_peptide_count + entrapment_peptide_count):.2f}%")
    print(f"sample method: {entrapment_peptide_count * (1 / r) * 100 / target_peptide_count:.2f}%")


def calculate(wd, entrapment_prefix, tool_name):
    df = pd.DataFrame()

    if tool_name == "fp":
        for p in glob.glob(os.path.join(wd, "**/peptide.tsv"), recursive=True):
            df2 = pd.read_csv(p, sep='\t', usecols=['Peptide', 'Protein', 'Mapped Proteins'], index_col=False)
            df2['all_proteins'] = df2.apply(lambda x: ','.join(filter(pd.notna, [x['Protein'], x['Mapped Proteins']])), axis=1)
            df2['IsEntrapment'] = df2['all_proteins'].apply(lambda pg: all(entrapment_prefix in p for p in pg.split(',')))
            df = pd.concat([df, df2], ignore_index=True)
    elif tool_name == "mm":
        for p in glob.glob(os.path.join(wd, "**-calib_Peptides.psmtsv"), recursive=True):
            df2 = pd.read_csv(p, sep="\t", index_col=False, na_values=["", "n/a"], header=0, usecols=["Base Sequence", "Contaminant", "Decoy", "PEP_QValue", "Protein Accession"])
            df2 = df2[df2["Contaminant"] == "N"]
            df2 = df2[df2["Decoy"] == "N"]
            df2 = df2[df2["PEP_QValue"] < 0.01]
            df2['IsEntrapment'] = df2['Protein Accession'].apply(lambda pg: all(entrapment_prefix in p for p in pg.split('|')))
            df2['Peptide'] = df2['Base Sequence']
            df = pd.concat([df, df2], ignore_index=True)
    elif tool_name == "mq":
        df = pd.read_csv(wd + r"/evidence.txt", low_memory=False, sep="\t", index_col=False, na_values=["", "NaN"], header=0, usecols=["Sequence", "Proteins", "Reverse", "Potential contaminant"])
        df = df[df["Reverse"].isna()]
        df = df[df["Potential contaminant"].isna()]
        df['IsEntrapment'] = df['Proteins'].apply(lambda pg: False if pd.isna(pg) else all(entrapment_prefix in p for p in pg.split(';')))
        df['Peptide'] = df['Sequence']
    else:
        print(f"Tool {tool_name} is not supported.")
        exit(1)

    target_peptide_count = df[~df['IsEntrapment']]['Peptide'].nunique(dropna=True)
    entrapment_peptide_count = df[df['IsEntrapment']]['Peptide'].nunique(dropna=True)

    return target_peptide_count, entrapment_peptide_count


if __name__ == "__main__":
    root = r'Z:/yufe/results/msfragger_ddaplus_paper/MSV000090552/'

    print("\nMSFragger-DDA+")
    main(["entrapment_", root + 'fragpipe_ddaplus_entrapment', 'fp', 1])

    print("\nMSFragger-DDA")
    main(["entrapment_", root + 'fragpipe_dda_entrapment', 'fp', 1])

    print("\nMetaMorpheus")
    # main(["entrapment_", root + r'metamorpheus_entrapment/Task2-SearchTask/Individual File Results/', 'mm', 1])
    main(["entrapment_", root + r'metamorpheus_entrapment_106/Task2-SearchTask/Individual File Results/', 'mm', 1])

    print("\nMaxQuant")
    main(["entrapment_", root + r'maxquant_entrapment/combined/txt', 'mq', 1])

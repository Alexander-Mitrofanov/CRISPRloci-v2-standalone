import os
import pandas as pd


def export_defense_finder_systems(defense_finder_genes, outdir, filename):
    systems = build_defense_finder_systems(defense_finder_genes)
    systems.to_csv(os.path.join(outdir, filename + '_defense_finder_systems.tsv'), sep='\t', index=False)


def build_defense_finder_systems(defense_finder_genes):
    """Summarize each system without a many-to-many merge on shared proteins.

    Modified 2026-10-12: a protein participating in multiple systems previously
    duplicated system rows during the final coordinate-ordering merge.
    """
    fields = ['sys_id', 'type', 'subtype', 'activity', 'sys_beg', 'sys_end',
              'protein_in_syst', 'genes_count', 'name_of_profiles_in_sys']
    rows = []
    for identity, group in defense_finder_genes.groupby(['sys_id', 'type', 'subtype', 'activity'], sort=True):
        ordered = group.sort_values(['hit_pos', 'hit_id'], kind='stable')
        proteins = list(dict.fromkeys(ordered.hit_id))
        row = dict(zip(fields[:4], identity))
        row.update(sys_beg=proteins[0], sys_end=proteins[-1], protein_in_syst=','.join(proteins),
                   genes_count=len(proteins), name_of_profiles_in_sys=','.join(sorted(group.gene_name)),
                   _position=int(ordered.hit_pos.iloc[0]))
        rows.append(row)
    rows.sort(key=lambda row: (row['_position'], row['sys_id']))
    return pd.DataFrame(rows, columns=fields)

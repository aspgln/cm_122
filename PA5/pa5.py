import numpy as np
import pandas as pd


if __name__ == "__main__":
    input_folder = 'hw5_E_1'
    input_fn_start = '{0}/{0}'.format(input_folder)
    exons_fn = '{}_exons.txt'.format(input_fn_start)
    isoforms_fn = '{}_isoforms.txt'.format(input_fn_start)
    exon_counts_fn = '{}_exon_counts.txt'.format(input_fn_start)

    exon_length_dict = {}
    with open(exons_fn) as exons_file:
        exons_file.readline()
        for line in exons_file:
            exon, start_end = line.strip().split(':')
            start, end = (int(x) for x in start_end.split(','))
            exon_length = end - start + 1
            exon_length_dict[exon] = exon_length
    # print exon_length_dict
    exon_length_df = pd.DataFrame.from_dict(exon_length_dict, orient='index').reset_index().rename(columns={'index':'exon', 0:'length'})
    print exon_length_df.head()

    isoform_exons = []
    with open(isoforms_fn) as isoforms_file:
        isoforms_file.readline()
        for line in isoforms_file:
            isoform, exons = line.strip().split(':')
            exons = exons.split(',')
            gene_id = isoform.split('_ISO')[0]
            full_exons = ['{}_{}'.format(gene_id, exon) for exon in exons]
            these_exons = [(isoform, exon) for exon in full_exons]
            isoform_exons += these_exons
    print isoform_exons

    isoforms_df = pd.DataFrame(isoform_exons).rename(columns={0:'isoform', 1:'exon'})
    print isoforms_df.head()

    isoform_exon_df = pd.merge(isoforms_df, exon_length_df, on='exon')
    print isoform_exon_df.head()

    combined_df = pd.pivot_table(isoform_exon_df, columns='isoform', index='exon', values='length').fillna(0)
    print combined_df.head()

    isoform_columns = combined_df.columns
    combined_df.reset_index(inplace=True)
    print combined_df.head()
    print combined_df.columns

    exon_counts = []
    with open(exon_counts_fn) as exon_counts_file:
        exon_counts_file.readline()
        for line in exon_counts_file:
            exon, count = line.strip().split(':')
            exon_counts.append((exon, int(count)))
    exon_counts_df = pd.DataFrame(exon_counts).rename(columns={0:'exon', 1:'exon_count'})
    print exon_counts_df.head()
    complete_df = pd.merge(combined_df, exon_counts_df, on='exon')
    print complete_df.head()

    exon_counts_array = complete_df.exon_count.values
    exon_isoform_matrix = complete_df.loc[:, isoform_columns].values
    print exon_counts_array.shape
    print exon_isoform_matrix.shape

    implied_frequencies = np.linalg.lstsq(exon_isoform_matrix, exon_counts_array)[0]
    implied_frequencies /= implied_frequencies.sum()
    print implied_frequencies

    output_fn = '{}_output.txt'.format(input_fn_start)
    with open(output_fn, 'w') as output_file:
        output_file.write('>{}\n'.format(input_folder))
        for isoform, freq in sorted([_ for _ in zip(isoform_columns, implied_frequencies)],
                                    key=lambda x: (int(x[0].split('_')[1]), int(x[0].split('_')[3]))):
            output_file.write('{}:{}\n'.format(isoform, freq))

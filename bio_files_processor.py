def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None) -> str:
    '''
    Converts multi-line FASTA sequences into single-line format and save result.

    Args:
        input_fasta: Path to input FASTA file.
        output_fasta: Optional path to save converted FASTA.

    Returns:
        Path to output FASTA file.
    '''

    if not output_fasta:
        file_format = input_fasta.split('.')[-1]
        output_fasta = input_fasta.replace(f'.{file_format}', f'_oneline.{file_format}')

    with open(input_fasta, 'r') as input_fasta_file, open(output_fasta, 'w') as output_fasta_file:
        reading_start = True
        for row in input_fasta_file:
            if row.startswith('>'):
                if reading_start:
                    output_fasta_file.write(row)
                    reading_start = False
                else:
                    output_fasta_file.write(f'\n{row}')
            else:
                output_fasta_file.write(row.strip())
        output_fasta_file.write('\n')
    return output_fasta


def parse_blast_output(input_file: str, output_file: str = None) -> str:
    '''
    Parses a BLAST output file and keeps only the best match for each query sequence.

    Args:
        input_file: Path to the BLAST output file.
        output_file: Optional path to save filtered results (best hits only).

    Returns:
        Path to output table file.
    '''

    if not output_file:
        file_format = input_file.split('.')[-1]
        output_file = input_file.replace(f'.{file_format}', f'_best_matches.{file_format}')

    with open(input_file, 'r') as input_file_, open(output_file, 'w') as output_file_:
        proteins = []
        description = False
        for row in input_file_:
            if row.startswith('Description'):
                description = True
                second_column_start = row.find('Name')
            elif description:
                protein = row[0:second_column_start].strip()
                proteins.append(protein)
                description = False

        for prot in sorted(proteins):
            output_file_.write(f'{prot}\n')
    return output_file


def select_genes_from_gbk_to_fasta(input_gbk: str, genes: str | list, n_before: int = 1, n_after: int = 1, output_fasta: str = None) -> str:
    '''
    Extracts the nearest flanking genes of a target gene from a GenBank file and saves them in FASTA format.

    Args:
        input_gbk: Path to the input GenBank file.
        genes: Target gene name or list of gene names.
        n_before: Number of upstream (preceding) genes to include.
        n_after: Number of downstream (following) genes to include.
        output_fasta: Optional path to save the FASTA file.

    Returns:
        Path to the generated FASTA file.
    '''

    if isinstance(genes, str):
        genes = [genes]

    if not output_fasta:
        output_fasta = input_gbk.replace('.gbk', '_flanked_genes.fasta')

    with open(input_gbk, 'r') as input_file_, open(output_fasta, 'w') as output_file_:
        feature_buffer_size = n_before + n_after + 1
        feature_buffer = []

        gene_match_timers = []
        gene_matches = []

        features_start = False
        feature_data = {}
        for row in input_file_:
            if row.startswith('FEATURES'):
                second_column_start = row.find('Location/Qualifiers')
                features_start = True
            elif row.startswith('ORIGIN'):
                features_start = False
            elif features_start:
                first_column = row[0:second_column_start].strip()
                second_column = row[second_column_start:].strip()

                if first_column:
                    if feature_data:
                        if {'gene', 'translation'} <= feature_data.keys():  # Добавляем feature только с нужными тэгами
                            if len(feature_buffer) == feature_buffer_size:  # Поддержание буфера нужной размерности
                                del feature_buffer[0]
                            feature_buffer.append(feature_data)

                            gene = feature_data.get('gene')
                            # Если текущий ген есть в списке genes, значит буфер уже содержит данные для n_before.
                            # Теперь нужно набрать данных для n_after. Для этого мы для каждого подходящего гена
                            # создаем таймер равный кол-ву генов, которые нужно набрать - n_after + 1.
                            # +1 потому, что в буфер также включается текущий ген (целевой).
                            # Таким образом алгоритм пройдет еще n_after + 1 и соберет данные перед записью.
                            if gene in genes:
                                gene_match_timers.append(n_after + 1)  # Добавляем таймер для каждого мэтча
                                gene_matches.append(gene)

                            # Каждый раз, когда встречается подходящий feature (в данном случает содержащий gene и translation),
                            # от таймера для всех мэтчей отнимается 1.
                            # Когда таймер дойдет до 0 произойдет "сброс бомбы", это означает, что необходимое кол-во feature для покрытия n_after
                            # набрано и данные для текущего мэтча можно записывать. Оставшиеся мэтчи продолжат набирать свои n_after.
                            gene_match_timers = [x - 1 for x in gene_match_timers]
                            completed_timers = 0
                            for gene_match_timer, target_gene in zip(gene_match_timers, gene_matches):
                                if gene_match_timer == 0:
                                    for feature in feature_buffer:
                                        matched_gene = feature['gene']
                                        matched_translation = feature['translation']
                                        if target_gene != matched_gene:
                                            output_file_.write(f'>{matched_gene}_{target_gene}|target_gene={target_gene}\n')
                                            output_file_.write(f'{matched_translation}\n')
                                    completed_timers += 1

                            if completed_timers:
                                del gene_match_timers[0:completed_timers]
                                del gene_matches[0:completed_timers]

                        feature_data = {}

                    feature_data = {
                            'feature': first_column,
                            'location': second_column
                        }
                else:
                    if second_column.startswith('/'):
                        tag, data = second_column.split('=')
                        tag = tag.strip('/')
                        data = data.strip('"')
                        feature_data[tag] = data
                    else:
                        feature_data[tag] += data
    return output_fasta
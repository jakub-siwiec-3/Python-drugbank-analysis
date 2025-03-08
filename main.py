import xml.etree.ElementTree as ET
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import random


# Funkcja pomocnicza
def replace_spaces_with_nl(text):
    return text.replace(' ', '\n')


# Podpunkt 1: Ekstrakcja podstawowych informacji o leku
def extract_drug_info(root, namespace):
    data1 = []
    for drug in root.findall('ns:drug', namespace):
        drug_id = drug.find('ns:drugbank-id[@primary="true"]', namespace)
        name = drug.find('ns:name', namespace)
        drug_type = drug.attrib['type']
        description = drug.find('ns:description', namespace)
        state = drug.find('ns:state', namespace)
        indication = drug.find('ns:indication', namespace)
        mechanism_of_action = drug.find('ns:mechanism-of-action', namespace)

        food_interactions = drug.find('ns:food-interactions', namespace)
        food_interaction_list = []
        if food_interactions is not None:
            food_interaction_list = [
                fi.text for fi in food_interactions.findall('ns:food-interaction', namespace)
            ]

        data1.append({
            'Drug ID': drug_id.text,
            'Name': name.text,
            'Type': drug_type,
            'Description': description.text if description is not None else None,
            'State': state.text if state is not None else None,
            'Indication': indication.text if indication is not None else None,
            'Mechanism of action': mechanism_of_action.text if mechanism_of_action is not None else None,
            'Food interactions': food_interaction_list
        })
    return pd.DataFrame(data1)


# Podpunkt 2: Ekstrakcja synonimów
def extract_synonyms(root, namespace):
    data2 = []
    for drug in root.findall('ns:drug', namespace):
        drug_primary_id = drug.find('ns:drugbank-id[@primary="true"]', namespace)
        name = drug.find('ns:name', namespace)
        synonyms = drug.find('ns:synonyms', namespace)
        synonyms_list = []
        if synonyms is not None:
            synonyms_list = [s.text for s in synonyms.findall('ns:synonym', namespace)]
        data2.append({
            'drugbank_id': drug_primary_id.text,
            'name': name.text,
            'synonyms': synonyms_list
        })
    return pd.DataFrame(data2)


# Przykładowa funkcja rysująca graf synonimów (użycie z podpunktu 2)
def draw_synonym_graph(drugbank_id, df_synonyms):
    row = df_synonyms[df_synonyms['drugbank_id'] == drugbank_id]
    if row.empty:
        print(f"DrugBank ID {drugbank_id} not found.")
        return

    synonyms = row.iloc[0]['synonyms']
    nodes = [drugbank_id] + synonyms
    complete_graph = nx.complete_graph(nodes)
    options = {
        'node_size': 16000,
        'node_color': 'lightgreen',
        'edge_color': 'gray',
        'font_size': 7,
        'font_weight': 'bold',
        'with_labels': True
    }
    plt.figure(figsize=(10, 10))
    nx.draw_shell(complete_graph, **options)
    plt.show()


# Podpunkt 3: Ekstrakcja informacji o produktach
def extract_products(root, namespace):
    data3 = []
    for drug in root.findall('ns:drug', namespace):
        drug_id = drug.find('ns:drugbank-id[@primary="true"]', namespace)
        products = drug.find('ns:products', namespace)
        if products is not None:
            for product in products.findall('ns:product', namespace):
                product_name = product.find('ns:name', namespace)
                labeller = product.find('ns:labeller', namespace)
                dpd_id = product.find('ns:dpd-id', namespace)
                ndc_product_code = product.find('ns:ndc-product-code', namespace)
                dosage_form = product.find('ns:dosage-form', namespace)
                route = product.find('ns:route', namespace)
                strength = product.find('ns:strength', namespace)
                country = product.find('ns:country', namespace)
                source = product.find('ns:source', namespace)

                data3.append({
                    'Drug ID': drug_id.text,
                    'Product Name': product_name.text if product_name is not None else None,
                    'Labeller': labeller.text if labeller is not None else None,
                    'DPD ID': dpd_id.text if dpd_id is not None else None,
                    'NDC Product Code': ndc_product_code.text if ndc_product_code is not None else None,
                    'Dosage Form': dosage_form.text if dosage_form is not None else None,
                    'Route': route.text if route is not None else None,
                    'Strength': strength.text if strength is not None else None,
                    'Country': country.text if country is not None else None,
                    'Source': source.text if source is not None else None
                })
    return pd.DataFrame(data3)


# Podpunkt 4: Ekstrakcja informacji o ścieżkach (pathways)
def extract_pathways(root, namespace):
    data4 = []
    data4_v2 = []
    for drug in root.findall('ns:drug', namespace):
        if drug.tag == '{http://www.drugbank.ca}drug':
            drug_id = drug.find('ns:drugbank-id[@primary="true"]', namespace).text
            name = drug.find('ns:name', namespace)
            pathways = drug.find('ns:pathways', namespace)
            if pathways is not None and len(pathways.findall('ns:pathway', namespace)) > 0:
                for pathway in pathways.findall('ns:pathway', namespace):
                    pathway_name = pathway.find('ns:name', namespace)
                    drugs_list = []
                    drugs = pathway.find('ns:drugs', namespace)
                    for drug_in_pathway in drugs.findall('ns:drug', namespace):
                        new_drug = (drug_in_pathway.find('ns:drugbank-id', namespace).text,
                                    drug_in_pathway.find('ns:name', namespace).text)
                        drugs_list.append(new_drug)
                    data4.append({
                        'Pathway Name': pathway_name.text,
                        'Drug List': drugs_list
                    })
                    data4_v2.append({
                        'Drug ID': drug_id,
                        'Drug Name': name.text,
                        'Pathway Name': pathway_name.text,
                    })
            else:
                data4_v2.append({
                    'Drug ID': drug_id,
                    'Drug Name': name.text,
                    'Pathway Name': None
                })
    df_pathways = pd.DataFrame(data4)
    df_pathways_v2 = pd.DataFrame(data4_v2)
    return df_pathways, df_pathways_v2


# Podpunkt 5: Rysowanie grafu dwudzielnego (bipartite) dla ścieżek
def draw_pathway_bipartite_graph(data4):
    pathways = [replace_spaces_with_nl(item['Pathway Name']) for item in data4]
    B = nx.Graph()
    B.add_nodes_from(pathways, bipartite=0)
    for pathway_drugs in data4:
        for drug in pathway_drugs['Drug List']:
            B.add_node(drug[1], bipartite=1)
            B.add_edge(replace_spaces_with_nl(pathway_drugs['Pathway Name']), drug[1])
    pos = nx.bipartite_layout(B, nodes=pathways)
    plt.figure(figsize=(10, 20))
    nx.draw(B, pos, with_labels=True, node_size=8000,
            node_color=['lightblue' if node in pathways else 'lightgreen' for node in B.nodes()],
            font_size=11, alpha=0.8, width=1.0)
    plt.show()


# Podpunkt 6: Rysowanie histogramu występowania leków w różnych ścieżkach
def plot_drug_histogram(data4):
    drug_counts = {}
    for drug_list in [item['Drug List'] for item in data4]:
        for drug in drug_list:
            if drug in drug_counts:
                drug_counts[drug] += 1
            else:
                drug_counts[drug] = 1
    df_6 = pd.DataFrame(list(drug_counts.items()), columns=['Drug', 'Count'])
    df_6['Drug'] = df_6['Drug'].apply(lambda d: d[1])  # wyciągamy tylko nazwę leku z krotki
    colors = [[random.random() for _ in range(3)] for _ in range(len(df_6))]
    plt.figure(figsize=(10, 7))
    plt.bar(df_6['Drug'], df_6['Count'], color=colors)
    plt.title('Histogram of Drug Occurrences in Different Pathways')
    plt.xlabel('Drug name')
    plt.ylabel('Count of Pathways')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()


# Podpunkt 7: Ekstrakcja informacji o targetach
def extract_targets(root, namespace):
    data7 = []
    for drug in root.findall('ns:drug', namespace):
        if drug.tag == '{http://www.drugbank.ca}drug':
            drug_id = drug.find('ns:drugbank-id[@primary="true"]', namespace).text
            targets = drug.find('ns:targets', namespace)
            if targets is not None:
                for target in targets.findall('ns:target', namespace):
                    target_id = target.find('ns:id', namespace)
                    polypeptide = target.find('ns:polypeptide', namespace)
                    if polypeptide is not None:
                        src_id = polypeptide.attrib['id']
                        src_src = polypeptide.attrib['source']
                        name = polypeptide.find('ns:name', namespace)
                        gene_name = polypeptide.find('ns:gene-name', namespace)
                        cel_loc = polypeptide.find('ns:cellular-location', namespace)
                        chromosome_loc = polypeptide.find('ns:chromosome-location', namespace)
                        external_id = polypeptide.find('ns:external-identifiers', namespace)
                        found_id = None
                        if external_id is not None:
                            for x in external_id.findall('ns:external-identifier', namespace):
                                res_names = x.findall('ns:resource', namespace)
                                if res_names is not None:
                                    for res_name in res_names:
                                        if res_name.text == 'GenAtlas':
                                            found_id = x.find('ns:identifier', namespace)
                                            break
                        data7.append({
                            'Drug ID': drug_id,
                            'Target ID': target_id.text,
                            'External source': src_src,
                            'External ID': src_id,
                            'Polypeptide name': name.text,
                            'Gene name': gene_name.text if gene_name is not None else None,
                            'Cellular location': cel_loc.text if cel_loc is not None else None,
                            'Chromosome': chromosome_loc.text if chromosome_loc is not None else None,
                            'Gen Atlas ID': found_id.text if found_id is not None else None,
                        })
    return pd.DataFrame(data7)


# Podpunkt 8: Wykres kołowy dla rozmieszczenia targetów
def plot_targets_distribution(df_targets):
    plot8_data = df_targets.groupby(['Cellular location']).size()
    cmap = plt.get_cmap('tab20')
    colors = [cmap(i) for i in range(len(plot8_data))]
    plt.figure(figsize=(20, 20))
    plt.pie(plot8_data, autopct='%1.1f%%', pctdistance=1.1, labels=None, startangle=220, colors=colors)
    plt.legend(labels=plot8_data.index, loc='lower left', bbox_to_anchor=(-0.1, -0.1))
    plt.title('Percentage Distribution of Targets in Cellular Locations', fontweight='bold', fontsize=26)
    plt.show()


# Podpunkt 9: Ekstrakcja grup leków i wykresy kołowe
def extract_drug_groups(root, namespace):
    data9 = []
    for drug in root.findall('ns:drug', namespace):
        if drug.tag == '{http://www.drugbank.ca}drug':
            drug_id = drug.find('ns:drugbank-id[@primary="true"]', namespace)
            groups = drug.find('ns:groups', namespace)
            group_list = []
            if groups is not None:
                group_list = [s.text for s in groups.findall('ns:group', namespace)]
            data9.append({
                'Drug ID': drug_id.text if drug_id is not None else None,
                'Approved': 'approved' in group_list,
                'Vet Approved': 'vet_approved' in group_list,
                'Withdrawn': 'withdrawn' in group_list,
                'Investigational': 'investigational' in group_list or 'experimental' in group_list,
            })
    return pd.DataFrame(data9)


def plot_drug_groups(df_9):
    app_and_not_with = df_9[(df_9['Approved'] == True) & (df_9['Withdrawn'] == False)]
    print(f'Approved and not withdrawn amount: {app_and_not_with.shape[0]}')
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))

    def tf_pie(df, column_name, ax):
        value_counts = df[column_name].value_counts()
        colors_tf = ['#90EE90' if value else '#FFCCCB' for value in value_counts.index]
        ax.pie(value_counts, labels=value_counts.index, autopct='%1.1f%%', startangle=90, colors=colors_tf)
        ax.set_title(f'{column_name}')

    tf_pie(df_9, 'Approved', axes[0, 0])
    tf_pie(df_9, 'Vet Approved', axes[0, 1])
    tf_pie(df_9, 'Withdrawn', axes[1, 0])
    tf_pie(df_9, 'Investigational', axes[1, 1])

    plt.tight_layout()
    plt.show()


# Podpunkt 10: Ekstrakcja szczegółowych informacji o targetach oraz interakcji lekowych
def extract_targets_details(root, namespace):
    data10 = []
    for drug in root.findall('ns:drug', namespace):
        if drug.tag == '{http://www.drugbank.ca}drug':
            drug_id = drug.find('ns:drugbank-id[@primary="true"]', namespace).text
            targets = drug.find('ns:targets', namespace)
            for target in targets.findall('ns:target', namespace):
                target_id = target.find('ns:id', namespace)
                target_name = target.find('ns:name', namespace)
                target_organism = target.find('ns:organism', namespace)
                actions_list = []
                actions = target.find('ns:actions', namespace)
                for action in actions.findall('ns:action', namespace):
                    actions_list.append(action.text)
                polypeptide = target.find('ns:polypeptide', namespace)
                if polypeptide is not None:
                    src_id = polypeptide.attrib['id']
                    name = polypeptide.find('ns:name', namespace)
                    general_function = polypeptide.find('ns:general-function', namespace)
                    specific_function = polypeptide.find('ns:specific-function', namespace)
                    go_classifiers = polypeptide.find('ns:go-classifiers', namespace)
                    go_classifiers_list = []
                    for go_classifier in go_classifiers.findall('ns:go-classifier', namespace):
                        go_classifiers_list.append((
                            go_classifier.find('ns:category', namespace).text,
                            go_classifier.find('ns:description', namespace).text
                        ))
                    data10.append({
                        'Drug ID': drug_id,
                        'Target ID': target_id.text,
                        'Target name': target_name.text,
                        'Target organism': target_organism.text,
                        'Actions': actions_list,
                        'Polypeptide': {
                            'Polypeptide ID': src_id,
                            'Polypeptide name': name.text,
                            'Polypeptide general function': general_function.text if general_function is not None else None,
                            'Polypeptide specific function': specific_function.text if specific_function is not None else None,
                            'Polypeptide go classifiers': go_classifiers_list
                        }
                    })
                else:
                    data10.append({
                        'Drug ID': drug_id,
                        'Target ID': target_id.text,
                        'Target name': target_name.text,
                        'Target organism': target_organism.text,
                        'Actions': actions_list,
                        'Polypeptide': None
                    })
    return pd.DataFrame(data10)


def extract_drug_interactions(root, namespace):
    data10_v2 = []
    for drug in root.findall('ns:drug', namespace):
        if drug.tag == '{http://www.drugbank.ca}drug':
            drug_id = drug.find('ns:drugbank-id[@primary="true"]', namespace).text
            name = drug.find('ns:name', namespace).text
            drug_interactions = drug.find('ns:drug-interactions', namespace)
            for interaction in drug_interactions.findall('ns:drug-interaction', namespace):
                sec_id = interaction.find('ns:drugbank-id', namespace)
                sec_name = interaction.find('ns:name', namespace)
                sec_description = interaction.find('ns:description', namespace)
                data10_v2.append({
                    'Drug 1 ID': drug_id,
                    'Drug 1 name': name,
                    'Drug 2 ID': sec_id.text,
                    'Drug 2 name': sec_name.text,
                    'Description': sec_description.text if sec_description is not None else None,
                })
    df_10_v2 = pd.DataFrame(data10_v2)
    return df_10_v2.drop_duplicates(subset=['Description'])


# Podpunkt 11: Praca z konkretnym genem
def build_gene_graph(specified_gene, df_targets, df_products):
    # Wyszukiwanie leków, które mają podany gen (w kolumnie 'Gene name' z df_targets)
    drugs_with_specified_gene = [item["Drug ID"] for item in df_targets.to_dict(orient='records') if
                                 item["Gene name"] == specified_gene]

    GG = nx.Graph()
    GG.add_node(specified_gene)
    for drug in drugs_with_specified_gene:
        GG.add_node(drug)
        GG.add_edge(specified_gene, drug)
    for drug in drugs_with_specified_gene:
        for product in df_products.to_dict(orient='records'):
            if product['Drug ID'] == drug:
                product_name = product['Product Name']
                if product_name is not None:
                    product_name = replace_spaces_with_nl(product_name)
                    GG.add_node(product_name)
                    GG.add_edge(drug, product_name)
    options = {
        'node_size': 13000,
        'node_color': 'lightgreen',
        'edge_color': 'gray',
        'font_size': 8,
        'font_weight': 'bold',
        'with_labels': True
    }
    plt.figure(figsize=(35, 35))
    nx.draw(GG, **options)
    plt.show()


def extract_gene_info(specified_gene, root, namespace):
    gene_info = {}
    gene_synonyms_list = []
    gene_ATGC = None
    for drug in root.findall('ns:drug', namespace):
        if drug.tag == '{http://www.drugbank.ca}drug':
            targets = drug.find('ns:targets', namespace)
            if targets is not None:
                for target in targets.findall('ns:target', namespace):
                    polypeptide = target.find('ns:polypeptide', namespace)
                    if polypeptide is not None:
                        gene_name = polypeptide.find('ns:gene-name', namespace)
                        if gene_name is not None and gene_name.text == specified_gene:
                            locus = polypeptide.find("ns:locus", namespace)
                            cellular_location = polypeptide.find("ns:cellular-location", namespace)
                            transmembrane_regions = polypeptide.find("ns:transmembrane-regions", namespace)
                            signal_regions = polypeptide.find("ns:signal-regions", namespace)
                            theoretical_pi = polypeptide.find("ns:theoretical-pi", namespace)
                            molecular_weight = polypeptide.find("ns:molecular-weight", namespace)
                            chromosome_location = polypeptide.find("ns:chromosome-location", namespace)
                            organism = polypeptide.find("ns:organism", namespace)
                            ncbi_taxonomy_id = organism.attrib["ncbi-taxonomy-id"] if organism is not None else None
                            gene_synonyms = polypeptide.find('ns:synonyms', namespace)
                            if gene_synonyms is not None:
                                for synonym in gene_synonyms.findall('ns:synonym', namespace):
                                    gene_synonyms_list.append(synonym.text)
                            gene_ATGC = polypeptide.find('ns:gene-sequence', namespace).text
                            gene_info = {
                                "Name": gene_name.text,
                                "Locus": locus.text if locus is not None else None,
                                "Cellular Location": cellular_location.text if cellular_location is not None else None,
                                "Transmembrane Regions": transmembrane_regions.text if transmembrane_regions is not None else None,
                                "Signal Regions": signal_regions.text if signal_regions is not None else None,
                                "Theoretical pi": float(theoretical_pi.text) if theoretical_pi is not None else None,
                                "Molecular Weight": float(
                                    molecular_weight.text) if molecular_weight is not None else None,
                                "Chromosome Location": chromosome_location.text if chromosome_location is not None else None,
                                "Organism": organism.text if organism is not None else None,
                                "NCBI Taxonomy ID": ncbi_taxonomy_id,
                                "Gene Synonyms": gene_synonyms_list
                            }
                            break
    return gene_info, gene_ATGC


def plot_gene_nucleotide_structure(gene_ATGC, specified_gene):
    # Zakładamy, że gene_ATGC zawiera nagłówek i sekwencję nukleotydów
    split_sequence = gene_ATGC.split("\n", 1)[1]
    no_whitespace_sequence = "".join(split_sequence.split())
    color_map = {'A': 'green', 'T': 'red', 'C': 'blue', 'G': 'yellow'}
    plt.figure(figsize=(10, 2))
    for i, base in enumerate(no_whitespace_sequence):
        plt.gca().add_patch(mpatches.Rectangle((i, 0), 1, 1, color=color_map.get(base, 'gray')))
    plt.xlim(0, len(no_whitespace_sequence))
    plt.ylim(0, 1)
    plt.yticks([])
    plt.xticks([])
    legend_labels = [mpatches.Patch(color=color_map[letter], label=letter) for letter in 'ATCG']
    plt.legend(handles=legend_labels, bbox_to_anchor=(1, 1))
    plt.title(f"DNA Sequence Representation of {specified_gene}", fontsize=14)
    plt.show()


# Podpunkt 12: Ekstrakcja klasyfikacji i wykresy słupkowe
def extract_classification(root, namespace):
    data12 = []
    for drug in root.findall('ns:drug', namespace):
        if drug.tag == '{http://www.drugbank.ca}drug':
            drug_id = drug.find('ns:drugbank-id[@primary="true"]', namespace).text
            classification = drug.find('ns:classification', namespace)
            if classification is not None:
                direct_parent = classification.find('ns:direct-parent', namespace)
                superclass = classification.find('ns:superclass', namespace)
                data12.append({
                    'Drug ID': drug_id,
                    'Direct Parent': direct_parent.text if direct_parent is not None else None,
                    'Superclass': superclass.text if superclass is not None else None,
                })
    return pd.DataFrame(data12)


def plot_classification(df12):
    direct_parent_counts = df12['Direct Parent'].value_counts()
    superclass_counts = df12['Superclass'].value_counts()

    plt.figure(figsize=(15, 5))

    plt.subplot(1, 2, 1)
    direct_parent_counts.plot(kind='bar', color='blue')
    plt.title('Number of Drugs by Direct Parent')
    plt.xlabel('Direct Parent')
    plt.ylabel('Count')
    plt.xticks(rotation=45, ha='right')

    plt.subplot(1, 2, 2)
    superclass_counts.plot(kind='bar', color='purple')
    plt.title('Number of Drugs by Superclass')
    plt.xlabel('Superclass')
    plt.ylabel('Count')
    plt.xticks(rotation=45, ha='right')

    plt.tight_layout()
    plt.show()


# Główna funkcja, która wywołuje wszystkie kroki
def main():
    tree = ET.parse('drugbank_partial.xml')
    root = tree.getroot()
    namespace = {'ns': 'http://www.drugbank.ca'}

    # Podpunkt 1
    df_inf = extract_drug_info(root, namespace)
    # print(df_inf.to_string())

    # Podpunkt 2
    df_synonyms = extract_synonyms(root, namespace)
    # print(df_synonyms.to_string())
    draw_synonym_graph('DB00001', df_synonyms)

    # Podpunkt 3
    df_products = extract_products(root, namespace)
    # print(df_products.to_string())

    # Podpunkt 4
    df_pathways, df_pathways_v2 = extract_pathways(root, namespace)
    print(f"Number of unique pathways: {df_pathways['Pathway Name'].nunique()}")

    # Podpunkt 5
    draw_pathway_bipartite_graph(df_pathways.to_dict(orient='records'))

    # Podpunkt 6
    plot_drug_histogram(df_pathways.to_dict(orient='records'))

    # Podpunkt 7
    df_targets = extract_targets(root, namespace)
    # print(df_targets.to_string())

    # Podpunkt 8
    plot_targets_distribution(df_targets)

    # Podpunkt 9
    df_9 = extract_drug_groups(root, namespace)
    plot_drug_groups(df_9)

    # Podpunkt 10
    df_targets_details = extract_targets_details(root, namespace)
    df_drug_interactions = extract_drug_interactions(root, namespace)
    # Możesz tutaj wykorzystać lub zapisać te tabele

    # Podpunkt 11
    specified_gene = 'IFNAR2'
    build_gene_graph(specified_gene, df_targets, df_products)
    gene_info, gene_ATGC = extract_gene_info(specified_gene, root, namespace)
    print(gene_info)
    plot_gene_nucleotide_structure(gene_ATGC, specified_gene)

    # Podpunkt 12
    df12 = extract_classification(root, namespace)
    plot_classification(df12)


if __name__ == "__main__":
    main()

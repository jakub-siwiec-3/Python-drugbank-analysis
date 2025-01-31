import xml.etree.ElementTree as ET

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import random as random


tree = ET.parse('drugbank_partial.xml')
root = tree.getroot()

namespace = {'ns': 'http://www.drugbank.ca'}

def replace_spaces_with_nl(text):
    return text.replace(' ', '\n')

data1 = []
# Znajdowanie tagów <drugbank-id primary="true">
# Przechodzimy przez wszystkie dzieci <drug> w <drugbank> i zapisujemy interesujące nas dane

#podpunkt 1
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
        'Food interactions': food_interaction_list if food_interaction_list is not None else None
    })

df_inf = pd.DataFrame(data1)

# print(df_inf.to_string())

# podpunkt 2
# zapisywanie synonimów
data2 = []
for drug in root.findall('ns:drug', namespace):
    drug_primary_id = drug.find('ns:drugbank-id[@primary="true"]', namespace)
    name = drug.find('ns:name', namespace)
    synonyms = drug.find('ns:synonyms', namespace)
    synonyms_list = []
    if synonyms is not None:
        synonyms_list = [
            s.text for s in synonyms.findall('ns:synonym', namespace)
        ]

    data2.append({
        'drugbank_id': drug_primary_id.text,
        'name': name.text,
        'synonyms': synonyms_list
    })
df_synonyms = pd.DataFrame(data2)

# print(df_synonyms.to_string())

# tworzenie grafu pełnego - krawędź gdy synonimy
def draw_synonym_graph(drugbank_id):
    row = df_synonyms[df_synonyms['drugbank_id'] == drugbank_id]
    if row.empty:
        print(f"DrugBank ID {drugbank_id} not found.")
        return

    synonyms = row.iloc[0]['synonyms']

    # wybrany sposób reprezentacji leku
    nodes = [drugbank_id] + synonyms
    # nodes = [name.text] + synonyms

    complete_graph = nx.complete_graph(nodes)
    options = {
        'node_size': 16000,
        'node_color': 'lightgreen',
        'edge_color': 'gray',
        'font_size': 7,
        'font_weight': 'bold',
        'with_labels': True
    }
    plt.figure(figsize=(10,10))

    nx.draw_shell(complete_graph, **options)

    plt.show()

# przykladowy graf
draw_synonym_graph('DB00001')

# podpunkt 3
data3 = []
# wyszukujemy odpowiednich sekcji w strukturze xml i zapisujemy dane o produktach
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
                'Product Name': product_name.text,
                'Labeller': labeller.text if labeller is not None else None,
                'DPD ID': dpd_id.text if dpd_id is not None else None,
                'NDC Product Code': ndc_product_code.text if ndc_product_code is not None else None,
                'Dosage Form': dosage_form.text if dosage_form is not None else None,
                'Route': route.text if route is not None else None,
                'Strength': strength.text if strength is not None else None,
                'Country': country.text if country is not None else None,
                'Source': source.text if source is not None else None
            })

df_products = pd.DataFrame(data3)
# przykładowe wypisania
# print(df_products.to_string())
#print(df_products[df_products['Drug ID'] == 'DB00001'].to_string())


# podpunkt 4
# pierwotnie zakładałem że chodzi o leki w tagu 'drugs' wewnątrz pathway
# może jednak chodzić także o lek, który ma w swoim poddrzewie dane 'pathways'
# dlatego są 2 wersje; graf przygotowałem tylko dla 1. opcji bo dla 2. nie ma kompletnie sensu
data4 = []
data4_v2 = []
for drug in root.findall('ns:drug', namespace):
    if drug.tag == '{http://www.drugbank.ca}drug':
        drug_id = drug.find('ns:drugbank-id[@primary="true"]', namespace)
        drug_id = drug_id.text
        name = drug.find('ns:name', namespace)

        pathways = drug.find('ns:pathways', namespace)
        v2_list = []
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
                    'Drug ID' : drug_id,
                    'Drug Name' : name.text,
                    'Pathway Name' : pathway_name.text,
                })
        else :
            data4_v2.append({
                'Drug ID' : drug_id,
                'Drug Name': name.text,
                'Pathway Name' : None
            })

df_pathways = pd.DataFrame(data4)
df_pathways_v2 = pd.DataFrame(data4_v2)
unique_pathways_count = df_pathways['Pathway Name'].nunique()
print(f"Number of unique pathways: {unique_pathways_count}")

# podpunkt 5

B = nx.Graph()
pathways = [replace_spaces_with_nl(pathway['Pathway Name']) for pathway in data4]
# graf dwudzielny, nazwy szlaków z lewej strony, leki z prawej
B.add_nodes_from(pathways, bipartite=0)
for pathway_drugs in data4:
    for drug in pathway_drugs['Drug List']:
        # drugi element krotki - nazwa
        B.add_node(drug[1], bipartite=1)
        #tworzymy krawędzie między lekami a ścieżkami, z którymi wchodzą w interakcje
        B.add_edge(replace_spaces_with_nl(pathway_drugs['Pathway Name']), drug[1])
pos = nx.bipartite_layout(B, nodes=pathways)
plt.figure(figsize=(10, 20))
#dostosowujemy rozmiar, kolory, obecność etykiet, przezroczystość (alfa)
nx.draw(
    B, pos, with_labels=True, node_size=8000,
    node_color=['lightblue' if node in pathways else 'lightgreen' for node in B.nodes()],
    font_size=11, alpha=0.8, width=1.0
)
plt.show()

# podpunkt 6
# 1. wersja
drug_counts = {}

drug_lists = [drug_list['Drug List'] for drug_list in data4]
for drug_list in drug_lists:
    for drug in drug_list:
        if drug in drug_counts:
            drug_counts[drug] += 1
        else:
            drug_counts[drug] = 1

df_6 = pd.DataFrame(list(drug_counts.items()), columns=['Drug', 'Count'])
# print((df_6.to_string()))
df_6['Drug'] = df_6['Drug'].apply(lambda only_name: only_name[1])
def random_color():
    return [random.random() for _ in range(3)]

# używamy losowych kolorów
colors = [random_color() for _ in range(len(df_6))]
plt.figure(figsize=(10, 7))
# ustawiamy kolory do słupków
plt.bar(df_6['Drug'], df_6['Count'], color=colors)

plt.title('Histogram of Drug Occurrences in Different Pathways')
plt.xlabel('Drug name')
plt.ylabel('Count of Pathways')

plt.xticks(rotation=45)

plt.tight_layout()
plt.show()

# 2. wersja
drug_counts_v2 = {drug_id: 0 for drug_id in df_pathways_v2['Drug Name'].unique()}
# iterujemy się po wierszach i aktualizujemy licznik
for _, row in df_pathways_v2.iterrows():
    if row['Pathway Name'] is not None:
        drug_counts_v2[row['Drug Name']] += 1
df_counts = pd.DataFrame(list(drug_counts_v2.items()), columns=['Drug Name', 'Pathway Count'])
plt.figure(figsize=(16, 8))
plt.bar(df_counts['Drug Name'], df_counts['Pathway Count'])

plt.title('Count of pathways for each Drug', fontsize=16)
plt.xlabel('Drug Name', fontsize=14)
plt.ylabel('Count of Pathways', fontsize=14)

plt.xticks(rotation=90)

plt.tight_layout()
plt.show()

# podpunkt 7
data7 = []

for drug in root.findall('ns:drug', namespace):
    if drug.tag == '{http://www.drugbank.ca}drug':
        drug_id = drug.find('ns:drugbank-id[@primary="true"]', namespace)
        drug_id = drug_id.text

        targets = drug.find('ns:targets', namespace)
        if targets is not None:
            for target in targets.findall('ns:target', namespace):
                target_id = target.find('ns:id', namespace)
                polypeptide = target.find('ns:polypeptide', namespace)
                if polypeptide is not None:
                    src_id=polypeptide.attrib['id']
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
                                    if res_name.text == 'GenAtlas' :
                                        found_id = x.find('ns:identifier', namespace)
                                        break
                                    else :
                                        continue
                    data7.append({
                        'Drug ID' : drug_id,
                        'Target ID' : target_id.text,
                        'External source' : src_src,
                        'External ID' : src_id,
                        'Polypeptide name' : name.text,
                        'Gene name' : gene_name.text if gene_name is not None else None,
                        'Cellular location' : cel_loc.text if chromosome_loc is not None else None,
                        'Chromosome' : chromosome_loc.text if chromosome_loc is not None else None,
                        'Gen Atlas ID' : found_id.text if found_id is not None else None,

                })


df_targets = pd.DataFrame(data7)
# print(df_targets.to_string())

# podpunkt 8
plot8_data = df_targets.groupby(['Cellular location']).size()

# pobieramy mapę 20 kolorów z plt
cmap = plt.get_cmap('tab20')
# przypisujemy je do targetów
colors = [cmap(i) for i in range(len(plot8_data))]
x = plt.subplots(figsize=(20,20))
# tworzymy wykres kołowy z legendą; chcemy same podpisy procentów, bo inaczej będzie nieczytelne
plt.pie(plot8_data, autopct='%1.1f%%', pctdistance=1.1, labels=None, startangle=220, colors=colors)
# ustawiamy przypisy i pozycję na lewy dół dodatkowo przesunięty jeszcze bardziej w lewo i dół
plt.legend(labels=plot8_data.index, loc='lower left', bbox_to_anchor=(-0.1, -0.1))

plt.title('Percentage Distribution of Targets in Cellular Locations', fontweight='bold', fontsize=26)

plt.show()


# podpunkt 9
data9 = []
# etykiety z groups w ramce
for drug in root.findall('ns:drug', namespace):
    if drug.tag == '{http://www.drugbank.ca}drug':
        drug_id = drug.find('ns:drugbank-id[@primary="true"]', namespace)
        drug_id = drug_id.text if drug_id is not None else None
        groups = drug.find('ns:groups', namespace)
        group_list=[]
        if groups is not None:
            group_list = [
                s.text for s in groups.findall('ns:group', namespace)
            ]
        data9.append({
            'Drug ID':drug_id,
            'Approved' : 'approved' in group_list,
            'Vet Approved' : 'vet_approved' in group_list,
            'Withdrawn' : 'withdrawn' in group_list,
            'Investigational' : 'investigational' in group_list or 'experimental' in group_list,
        })

df_9 = pd.DataFrame(data9)
# print(df_9.to_string())

# liczba wierszy - wartosci spełniających warunek: zatwierdzone, nie wycofane
app_and_not_with = df_9[(df_9['Approved'] == True) & (df_9['Withdrawn'] == False)]
print(f'Approved and not withdrawn amount: {app_and_not_with.shape[0]}')

# wykresy kołowe dla każdego z kryteriów
def tf_pie(df, column_name, ax):
    # tworzymy ramkę danych z liczbą pozytywnych wartości dla każdego kryterium
    value_counts = df[column_name].value_counts()
    colors_tf = ['#90EE90' if value else '#FFCCCB' for value in value_counts.index]
    ax.pie(value_counts, labels=value_counts.index, autopct='%1.1f%%', startangle=90, colors=colors_tf)
    ax.set_title(f'{column_name}')

# tworzymy 4 podwykresy
fig, axes = plt.subplots(2, 2, figsize=(12, 12))

tf_pie(df_9, 'Approved', axes[0, 0])
tf_pie(df_9, 'Vet Approved', axes[0, 1])
tf_pie(df_9, 'Withdrawn', axes[1, 0])
tf_pie(df_9, 'Investigational', axes[1, 1])

plt.tight_layout()
plt.show()

# podpunkt 10
# zaszła w nim zmiana, więc zostawiam 2 wersje
data10 = []

for drug in root.findall('ns:drug', namespace):
    if drug.tag == '{http://www.drugbank.ca}drug':
        drug_id = drug.find('ns:drugbank-id[@primary="true"]', namespace)
        drug_id = drug_id.text
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
                    go_classifiers_list.append((go_classifier.find('ns:category', namespace).text,
                                                go_classifier.find('ns:description', namespace).text))
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
            else :
                data10.append({
                    'Drug ID': drug_id,
                    'Target ID': target_id.text,
                    'Target name': target_name.text,
                    'Target organism': target_organism.text,
                    'Actions': actions_list,
                    'Polypeptide' : None
                })

df_10 = pd.DataFrame(data10)
# print(df_10.to_string())

data10_v2 = []
# znajdujemy pary leków które ze sobą wchodzą w interakcje
# usuniemy tylko te krotki, dla których 'Description' są sobie równe
for drug in root.findall('ns:drug', namespace):
    if drug.tag == '{http://www.drugbank.ca}drug':
        drug_id = drug.find('ns:drugbank-id[@primary="true"]', namespace)
        drug_id = drug_id.text
        name = drug.find('ns:name', namespace).text
        drug_interactions = drug.find('ns:drug-interactions', namespace)
        for interaction in drug_interactions.findall('ns:drug-interaction', namespace):
            sec_id = interaction.find('ns:drugbank-id', namespace)
            sec_name = interaction.find('ns:name', namespace)
            sec_description = interaction.find('ns:description', namespace)
            data10_v2.append({
                'Drug 1 ID': drug_id,
                'Drug 1 name': name,
                'Drug 2 ID': sec_id,
                'Drug 2 name': sec_name,
                'Description' : sec_description.text if sec_description is not None else None,
            })
df_10_v2 = pd.DataFrame(data10_v2)
df_10_v2_filtered = df_10_v2.drop_duplicates(subset=['Description'])

# podpunkt 11 dla 1 konkretnego genu, zawartość:
# - graf prezentujący leki, które wchodzą w interakcję z tym genem oraz produkty które zawierają te leki
# - słownik z informajcami podstawowymi nt. genu
# - grafikę prezentującą strukturę nukleotydu

specified_gene = 'IFNAR2'
drugs_with_specified_gene = [item["Drug ID"] for item in data7 if item["Gene name"] == specified_gene]

# graf z krawędziami: gen-leki oraz leki-(produkty zawierające te leki)
GG = nx.Graph()
GG.add_node(specified_gene)
for drug in drugs_with_specified_gene:
    GG.add_node(drug)
    GG.add_edge(specified_gene, drug)
for drug in drugs_with_specified_gene:
    for drug_data in data3:
        if drug_data['Drug ID'] == drug:
            product_name = drug_data['Product Name']
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

# informacje o genie
gene_info = {}
gene_synonyms_list = []
# napis zawierający litery ATCG dla adeniny, guaniny, tyminy, cytozyny
gene_ATGC = None
for drug in root.findall('ns:drug', namespace):
    if drug.tag == '{http://www.drugbank.ca}drug':
        drug_id = drug.find('ns:drugbank-id[@primary="true"]', namespace)
        drug_id = drug_id.text

        targets = drug.find('ns:targets', namespace)
        if targets is not None:
            for target in targets.findall('ns:target', namespace):
                target_id = target.find('ns:id', namespace)
                polypeptide = target.find('ns:polypeptide', namespace)
                if polypeptide is not None:
                    gene_name = polypeptide.find('ns:gene-name', namespace)
                    if gene_name.text == specified_gene:
                        # pobieramy ogólne informacje o genie
                        locus = polypeptide.find("ns:locus", namespace)
                        cellular_location = polypeptide.find("ns:cellular-location", namespace)
                        transmembrane_regions = polypeptide.find("ns:transmembrane-regions", namespace)
                        signal_regions = polypeptide.find("ns:signal-regions", namespace)
                        theoretical_pi = polypeptide.find("ns:theoretical-pi", namespace)
                        molecular_weight = polypeptide.find("ns:molecular-weight", namespace)
                        chromosome_location = polypeptide.find("ns:chromosome-location", namespace)
                        organism = polypeptide.find("ns:organism", namespace)
                        ncbi_taxonomy_id=None
                        if organism is not None:
                            ncbi_taxonomy_id = organism.attrib["ncbi-taxonomy-id"]
                        gene_synonyms = polypeptide.find('ns:synonyms', namespace)
                        for synonym in gene_synonyms.findall('ns:synonym', namespace):
                            gene_synonyms_list.append(synonym.text)
                        gene_ATGC = polypeptide.find('ns:gene-sequence', namespace).text
                        gene_info = {
                            "Name" : gene_name.text,
                            "Locus": locus.text if locus is not None else None,
                            "Cellular Location": cellular_location.text if cellular_location is not None else None,
                            "Transmembrane Regions": transmembrane_regions.text if transmembrane_regions is not None else None,
                            "Signal Regions": signal_regions.text if signal_regions is not None else None,
                            "Theoretical pi": float(theoretical_pi.text) if theoretical_pi is not None else None,
                            "Molecular Weight": float(molecular_weight.text) if molecular_weight is not None else None,
                            "Chromosome Location": chromosome_location.text if chromosome_location is not None else None,
                            "Organism": organism.text if organism is not None else None,
                            "NCBI Taxonomy ID": ncbi_taxonomy_id if ncbi_taxonomy_id is not None else None,
                            "Gene Synonyms": gene_synonyms_list
                        }
                        break
                    else:
                        continue
# print(gene_info)


# tworzenie grafiki z nukleotydami

# ignorujemy pierwszy wiersz dzieląc napis raz i biorąc drugą część
split_sequence = gene_ATGC.split("\n", 1)[1]
# usuwamy białe spacje (domyślne split()) i mamy sekwencję nukleotydów
no_whitespace_sequence = "".join(split_sequence.split())

#tworzymy grafikę struktury nukleotydów
color_map = {'A': 'green', 'T': 'red', 'C': 'blue', 'G': 'yellow'}

plt.figure(figsize=(10, 2))

# rysujemy prostokąty dla nukleotydów
for i, base in enumerate(no_whitespace_sequence):
    plt.gca().add_patch(mpatches.Rectangle((i, 0), 1, 1, color=color_map[base]))

plt.xlim(0, len(no_whitespace_sequence))
plt.ylim(0, 1)

# usuwamy podpisy
plt.yticks([])
plt.xticks([])

# używamy mpatches do edycji każdego z elementów legendy - 4 obiekty Patch
legend_labels = [mpatches.Patch(color=color_map[letter], label=letter) for letter in 'ATCG']
plt.legend(handles=legend_labels, bbox_to_anchor=(1, 1))

plt.title(f"DNA Sequence Representation of {specified_gene}", fontsize=14)

plt.show()

#podpunkt 12
data12 = []
for drug in root.findall('ns:drug', namespace):
    if drug.tag == '{http://www.drugbank.ca}drug':
        drug_id = drug.find('ns:drugbank-id[@primary="true"]', namespace)
        drug_id = drug_id.text
        classification = drug.find('ns:classification', namespace)
        if classification is not None:
            direct_parent = classification.find('ns:direct-parent', namespace)
            kingdom = classification.find('ns:kingdom', namespace)
            superclass = classification.find('ns:superclass', namespace)
            data12.append({
                'Drug ID' : drug_id,
                'Direct Parent' : direct_parent.text if direct_parent is not None else None,
                'Kingdom' : kingdom.text if kingdom is not None else None,
                'Superclass' : superclass.text if superclass is not None else None,
            })
df12 = pd.DataFrame(data12)

# tworzymy wykresy przedstawiające liczebność różnych leków w poszczególnych klasyfikacjach
direct_parent_counts = df12['Direct Parent'].value_counts()
kingdom_counts = df12['Kingdom'].value_counts()
superclass_counts = df12['Superclass'].value_counts()

plt.figure(figsize=(15, 5))

# 1 wiersz, 3 kolumny i pierwszy indekx na tym wykresie
plt.subplot(1, 3, 1)
direct_parent_counts.plot(kind='bar', color='blue')
plt.title('Number of Drugs by Direct Parent')
plt.xlabel('Direct Parent')
plt.ylabel('Count')
plt.xticks(rotation=45, ha='right')

plt.subplot(1, 3, 2)
kingdom_counts.plot(kind='bar', color='green')
plt.title('Number of Drugs by Kingdom')
plt.xlabel('Kingdom')
plt.ylabel('Count')
plt.xticks(rotation=45, ha='right')

plt.subplot(1, 3, 3)
superclass_counts.plot(kind='bar', color='purple')
plt.title('Number of Drugs by Superclass')
plt.xlabel('Superclass')
plt.ylabel('Count')
plt.xticks(rotation=45, ha='right')

plt.tight_layout()
plt.show()
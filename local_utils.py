import numpy as np
import pandas as pd
import requests
from scipy import stats

def compute_mean(df, *args):
    """
    Function that takes 1 argument that is the dataframe and
    2 or 3 arguments that represent column names of given dataframe
    and computes the mean of each row for the given columns.
    """
    mean_list = []
    for row in df.iterrows():
        
        if len(args) == 3:
            col1 = args[0]
            col2 = args[1]
            col3 = args[2]
        
            #access values corresponding to each row for the given 3 columns
            values = df[[col1, col2, col3]].loc[row[0]].values
            zeros = len(values) - np.count_nonzero(values)

            if zeros == 0:      
                mean = values.mean()
                mean_list.append(mean)
            if zeros == 1:
                selected_values = values[values != 0] # find elements that are non-zero
                mean = selected_values.mean()
                mean_list.append(mean)

            if zeros == 2 or zeros == 3: #in case there is only one non-zero value or if all the values are zero
                mean = 0
                mean_list.append(mean)
        
        if len(args) == 2:
            col1 = args[0]
            col2 = args[1]
            
            values = df[[col1, col2]].loc[row[0]].values
            zeros = len(values) - np.count_nonzero(values)

            if zeros == 0:      
                mean = values.mean()
                mean_list.append(mean)

            if zeros == 1:
                selected_values = values[values != 0] # find elements that are non-zero
                mean = selected_values[0]
                mean_list.append(mean)

            if zeros == 2: #in case there is only one non-zero value or if all the values are zero
                mean = 0
                mean_list.append(mean)
            
    return mean_list

def get_mean(df):
    """
    Function that imports the datafame and separates the individual replicates within one sample
    into a list of lists (e.g. [[G7S_1,G7S_2,G7S_3],[...],[...],...]).
    It returns the dataframe with the mean of each replicate.
    """
    import math
    
    #because first column is the gene name that we don't need for parsing right now 
    list_col = list(df)[1::]
    
    #store elements in a list of lists
    ll_col = [list_col[3*i:3*i+3]  for i in range(0,math.ceil(len(list_col)/3)) ]

    l = []
    count = 4
    for elem in ll_col:
        mean = compute_mean(df, elem[0], elem[1], elem[2])
        df.insert(loc=count, column=str(elem[0][:-1])+'mean', value=mean)
        #elem.append(str(elem[0][:-1])+'mean')  
        count +=4

    return df
    
def get_mean2(df):
    """
    Computes meand of the gene set for sample Y3S that has 2 replicates instead of 3.
    """
    
    df["Y3S_mean"] = compute_mean(df, 'Y3S_1','Y3S_3')

    return df

def pca_comp(df, sorting='Cell state'):
    """
    Function that inputs a dataframe with genes as rows and features as columns 
    and outputs a PCA analysis with the genes given in the dataframe.
    """
    
    #dropping the first column that cointains genedID/symbol and storing it in another df
    df = df.iloc[:, 1:]
    
    #renaming the cells so that they all are labeled according to cell type or cell state
    if sorting == 'Cell state':
        df = df.rename(columns={"G3S": "Basal", "G7S": "Basal","G3F": "Basal", "G7F": "Basal","G3M": "Basal", "G7M": "Basal","Y3S": "Basal", "Y7S": "Basal","Y3F": "Basal", "Y7F": "Basal","Y3M": "Basal", "Y7M": "Basal", "G3SNSen": "NSen", "G7SNSen": "NSen","G3FNSen": "NSen", "G7FNSen": "NSen","G3MNSen": "NSen", "G7MNSen": "NSen","Y3SNSen": "NSen", "Y7SNSen": "NSen","Y3FNSen": "NSen", "Y7FNSen": "NSen","Y3MNSen": "NSen", "Y7MNSen": "NSen", "G3SSen": "Sen", "G7SSen": "Sen","G3FSen": "Sen", "G7FSen": "Sen","G3MSen": "Sen", "G7MSen": "Sen","Y3SSen": "Sen", "Y7SSen": "Sen","Y3FSen": "Sen", "Y7FSen": "Sen","Y3MSen": "Sen", "Y7MSen": "Sen"})
    else:
        df = df.rename(columns={"G3S": "Satellite", "G7S": "Satellite","G3F": "FAP", "G7F": "FAP","G3M": "Macrophage", "G7M": "Macrophage","Y3S": "Satellite", "Y7S": "Satellite","Y3F": "FAP", "Y7F": "FAP","Y3M": "Macrophage", "Y7M": "Macrophage", "G3SNSen": "Satellite", "G7SNSen": "Satellite","G3FNSen": "FAP", "G7FNSen": "FAP","G3MNSen": "Macrophage", "G7MNSen": "Macrophage","Y3SNSen": "Satellite", "Y7SNSen": "Satellite","Y3FNSen": "FAP", "Y7FNSen": "FAP","Y3MNSen": "Macrophage", "Y7MNSen": "Macrophage", "G3SSen": "Satellite", "G7SSen": "Satellite","G3FSen": "FAP", "G7FSen": "FAP","G3MSen": "Macrophage", "G7MSen": "Macrophage","Y3SSen": "Satellite", "Y7SSen": "Satellite","Y3FSen": "FAP", "Y7FSen": "FAP","Y3MSen": "Macrophage", "Y7MSen": "Macrophage"})

    #transposing the dataframe so  that the genes are columns and the conditions are rows
    genes_data = df.T
    
    # standardizing the data
    from sklearn.preprocessing import StandardScaler #transform data so distribution has mean = 0 and st.dev = 1
    scaler = StandardScaler()
    genes_data = scaler.fit_transform(genes_data)

    from sklearn.decomposition import PCA
    pca = PCA(n_components=2).fit(genes_data)
    pca_2d = pca.transform(genes_data)    

    #creating dataframe with umap values
    pca_df = pd.DataFrame(data = pca_2d , columns = ['PCA1', 'PCA2'])
    pca_df[sorting] = list(df.T.index)

    return(pca_df)

def kstest_function(series):
    """
    Function that applies to the dataframe that has sen first and nsen ans basal second in the oder
    (like the test df). It makes a KS test and appends the pval as last column to the dataframe.
    """
    series = series.copy()
    sen = series.iloc[0:36].to_numpy()
    nsen = series.iloc[36:].to_numpy()
    pval = stats.ks_2samp(sen,nsen)[1]
    return pval



def get_gene_ontology(gene_symbols, species="mouse"):
    """
    Queries MyGene.info to get Gene Ontology (GO) annotations.
    This is generally the simplest API for gene-centric data.
    """
    print(f"Querying Gene Ontology for: {', '.join(gene_symbols)}...")
    
    # MyGene.info API endpoint
    url = "https://mygene.info/v3/query"
    
    # Parameters for the query:
    # q: The query string (genes joined by comma)
    # scopes: Where to look for these strings (e.g., symbol, entrezgene)
    # fields: What data to return (we want 'go' for Gene Ontology)
    # species: 'mouse', 'human', 'rat', etc.
    params = {
        "q": ",".join(gene_symbols),
        "scopes": "symbol",
        "fields": "go",
        "species": species
    }
    
    try:
        response = requests.post(url, data=params)
        response.raise_for_status()
        data = response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data: {e}")
        return

    print("-" * 80)
    
    for item in data:
        symbol = item.get("query", "Unknown")
        go_data = item.get("go", {})
        
        # Check if GO data was found
        if not go_data:
            print(f"Gene: {symbol} - No GO annotations found.")
            print("-" * 80)
            continue

        print(f"Gene: {symbol}")
        
        # GO categories: 
        # BP = Biological Process
        # MF = Molecular Function
        # CC = Cellular Component
        categories = {"BP": "Biological Process", "MF": "Molecular Function", "CC": "Cellular Component"}
        
        for cat_code, cat_name in categories.items():
            terms = go_data.get(cat_code, [])
            
            # Sometimes terms is a dict (if only one result), normalize to list
            if isinstance(terms, dict):
                terms = [terms]
            if terms:
                print(f"  {cat_name}:")
                # OLD: Print top 3 terms as an example (to avoid spamming the console)
                termlist = []
                for term in terms:
                    if term.get('term', '') not in termlist:
                        term_id = term.get('id', '')
                        term_term = term.get('term', '')
                        print(f"    - [{term_id}] {term_term}")
                    termlist.append(term.get('term', ''))
                #if len(terms) > 3:
                #    print(f"    - ... and {len(terms) - 3} more")
        
        print("-" * 80)


def pca_variance(df):   
    #dropping the first column that cointains genedID/symbol and storing it in another df
    df = df.iloc[:, 1:]
    df = df.rename(columns={"G3S": "Basal", "G7S": "Basal","G3F": "Basal", "G7F": "Basal","G3M": "Basal", "G7M": "Basal","Y3S": "Basal", "Y7S": "Basal","Y3F": "Basal", "Y7F": "Basal","Y3M": "Basal", "Y7M": "Basal", "G3SNSen": "NSen", "G7SNSen": "NSen","G3FNSen": "NSen", "G7FNSen": "NSen","G3MNSen": "NSen", "G7MNSen": "NSen","Y3SNSen": "NSen", "Y7SNSen": "NSen","Y3FNSen": "NSen", "Y7FNSen": "NSen","Y3MNSen": "NSen", "Y7MNSen": "NSen", "G3SSen": "Sen", "G7SSen": "Sen","G3FSen": "Sen", "G7FSen": "Sen","G3MSen": "Sen", "G7MSen": "Sen","Y3SSen": "Sen", "Y7SSen": "Sen","Y3FSen": "Sen", "Y7FSen": "Sen","Y3MSen": "Sen", "Y7MSen": "Sen"})
    #transposing the dataframe so  that the genes are columns and features are rows
    genes_data = df.T
    
    # standardizing the data
    from sklearn.preprocessing import StandardScaler #transform data so distribution has mean = 0 and st.dev = 1
    scaler = StandardScaler()
    genes_data = scaler.fit_transform(genes_data)

    from sklearn.decomposition import PCA
    pca = PCA().fit(genes_data)
 
    # individual_variance holds the percentage of variance explained by each component
    ind_variance = pca.explained_variance_ratio_
    loadings = pca.components_

    return([ind_variance, loadings])

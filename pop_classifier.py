from sklearn.preprocessing import StandardScaler
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import cross_val_score
import sklearn.svm as svm
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse as ap
import os, sys

# Results in a warning that isn't useful in this context. Disabling.
pd.options.mode.chained_assignment = None  # default='warn'

def get_continental_mapping():
    """
    Returns a dictionary mapping population codes to continental groups.
    
    Returns:
        dict: A dictionary with population codes as keys and continental groups as values.
    """
    return {
        # East Asian (EAS)
        'CDX': 'EAS', 'CHB': 'EAS', 'CHS': 'EAS', 'JPT': 'EAS', 'KHV': 'EAS', 'CHD': 'EAS',
        
        # European (EUR)
        'CEU': 'EUR', 'FIN': 'EUR', 'GBR': 'EUR', 'IBS': 'EUR', 'TSI': 'EUR',
        
        # African (AFR)
        'ACB': 'AFR', 'ASW': 'AFR', 'ESN': 'AFR', 'GWD': 'AFR', 
        'LWK': 'AFR', 'MSL': 'AFR', 'YRI': 'AFR', 'MKK': 'AFR',
        
        # South Asian (SAS)
        'BEB': 'SAS', 'GIH': 'SAS', 'ITU': 'SAS', 'PJL': 'SAS', 'STU': 'SAS',
        
        # American (AMR)
        'CLM': 'AMR', 'MXL': 'AMR', 'PEL': 'AMR', 'PUR': 'AMR', 'MEX': 'AMR'
    }


def convert_to_continental(labels):
    """
    Convert population labels to continental groups.
    
    Args:
        labels (numpy.ndarray): Array of population labels.
        
    Returns:
        numpy.ndarray: Array of continental group labels.
    """
    continental_mapping = get_continental_mapping()
    
    # Find unique populations in the data
    unique_pops = np.unique(labels)
    
    # Check if there are any populations not in our mapping
    unknown_pops = [pop for pop in unique_pops if pop not in continental_mapping]
    if unknown_pops:
        print(f"Warning: Found populations not in continental mapping: {unknown_pops}")
        print("These will be kept as-is in the continental analysis.")
    
    # Convert to continental groups, keeping unknown populations as-is
    return np.array([continental_mapping.get(pop, pop) for pop in labels])


# First version that only returns superpopulation 

def run_old(pca_file, pop_file, linear=False):
    #parser = ap.ArgumentParser()
    #parser.add_argument('mergedpca', type=str, help='The merged pca file from the prePRIMUS pipeline.')
    #parser.add_argument('kgpops', type=str, help='The id/pop file for the 1000 genomes phase 3 samples')
    #parser.add_argument('--linear', action='store_true', help='Use LinearSVC(penalty="l2") instead of rbf #kernel. May result in results that make more sense but lower accuracy.')
    #args = parser.parse_args()

    # Load the 1000 Genomes PCA data
    all_pca = pd.read_csv(pca_file, header=None, sep=' ')
    all_pca.columns = ['FID', 'IID'] + [f'PC{i}' for i in range(1, 21)]

    kg_true = pd.read_csv(pop_file, sep="\t", header=None)
    kg_true.columns = ['IID', 'POP']

    # Convert 1000 Genomes true population labels to continental labels
    kg_true['continental'] = convert_to_continental(kg_true['POP'].values)

    # Extract rows from all_pca where IID is in kg_true IID column
    kg_data = all_pca[all_pca['IID'].isin(kg_true['IID'])]
    kg_data = kg_data.merge(kg_true, on='IID')

    # Extract PCs 1-10
    X_train = kg_data.iloc[:, 2:12].values
    y_train = kg_data['continental']

    # Probably not necessary because of the way the data is prepared but better safe than sorry
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)

    # From testing these were the best performing settings using HM3 as the train set and 1kg as the test.
    if linear_option == True:
        clf = svm.LinearSVC(penalty='l2')
    else:
        clf = svm.SVC(kernel='rbf', C=5, gamma=0.5)
    clf = CalibratedClassifierCV(clf, cv=5)  # Use cross-validation for calibration
    clf.fit(X_train_scaled, y_train)

    # Perform 5-fold cross-validation and calculate accuracy
    cv_scores = cross_val_score(clf, X_train_scaled, y_train, cv=5, scoring='accuracy')
    #print(f"5-Fold Cross-Validation Accuracy: {np.mean(cv_scores):.4f} ± {np.std(cv_scores):.4f}")

    # Get true continental and predicted continental categories for kg_data
    kg_data['predicted_continental'] = clf.predict(X_train_scaled)

    # Plot the first two principal components for 1000 Genomes data colored by population
    plt.figure(figsize=(10, 8))

    # Plot each population with a unique color
    populations = kg_data['continental'].unique()
    for pop in populations:
        subset = kg_data[kg_data['continental'] == pop]
        plt.scatter(subset['PC1'], subset['PC2'], label=pop, alpha=0.6)


    # Get directory to save figures 
    pca_file_bits = pca_file.split('/')[:-1]
    outpath = '/'.join(pca_file_bits)


    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('1000 Genomes PCA: Colored by Continental Population')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'{outpath}/compadre_1kg_pca_by_population.png')
    plt.close()

    # Plot the first two principal components for 1000 Genomes data
    plt.figure(figsize=(10, 8))

    # Plot all correctly classified samples in grey
    correct_samples = kg_data[kg_data['continental'] == kg_data['predicted_continental']]
    plt.scatter(correct_samples['PC1'], correct_samples['PC2'], label='Correctly Classified', color='grey', alpha=0.6)

    # Highlight misclassified samples with distinct colors
    misclassified_samples = kg_data[kg_data['continental'] != kg_data['predicted_continental']]
    for label in np.unique(misclassified_samples['predicted_continental']):
        subset = misclassified_samples[misclassified_samples['predicted_continental'] == label]
        plt.scatter(subset['PC1'], subset['PC2'], label=f'Misclassified as {label}', alpha=0.6, edgecolor='red')

    # Annotate misclassified samples with the true population
    for _, row in misclassified_samples.iterrows():
        plt.annotate(row['continental'], (row['PC1'], row['PC2']), fontsize=8, color='black')

    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('1000 Genomes PCA: True Labels and Misclassifications')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'{outpath}/compadre_1kg_pca_misclassifications.png')
    plt.close()

    # Plot higher PC comparisons (3:4, 5:6, etc.)
    plt.figure(figsize=(12, 10))

    pairs = [(3, 4), (5, 6), (7, 8), (9, 10)]

    for i, (pc_x, pc_y) in enumerate(pairs, start=1):
        plt.subplot(2, 2, i)

        # Plot each population with a unique color
        for pop in populations:
            subset = kg_data[kg_data['continental'] == pop]
            plt.scatter(subset[f'PC{pc_x}'], subset[f'PC{pc_y}'], label=pop, alpha=0.6)

        plt.xlabel(f'PC{pc_x}')
        plt.ylabel(f'PC{pc_y}')
        plt.grid(True)

    plt.suptitle('Higher PC Comparisons', fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Create a single legend for all subplots
    handles, labels = plt.gca().get_legend_handles_labels()
    unique_labels = dict(zip(labels, handles))
    plt.figlegend(unique_labels.values(), unique_labels.keys(), loc='upper center', ncol=3, fontsize='small', bbox_to_anchor=(0.5, 0.02))

    plt.savefig(f'{outpath}/compadre_1kg_pca_higher_pc.png', bbox_inches='tight')
    plt.close()

    # Plot higher PC comparisons (3:4, 5:6, etc.) with misclassification
    plt.figure(figsize=(12, 10))

    for i, (pc_x, pc_y) in enumerate(pairs, start=1):
        plt.subplot(2, 2, i)

        # Plot all correctly classified samples in grey
        plt.scatter(correct_samples[f'PC{pc_x}'], correct_samples[f'PC{pc_y}'], 
                    label='Correctly Classified', 
                    color='grey', alpha=0.6)

        # Highlight misclassified samples with distinct colors
        for label in np.unique(misclassified_samples['predicted_continental']):
            subset = misclassified_samples[misclassified_samples['predicted_continental'] == label]
            plt.scatter(subset[f'PC{pc_x}'], subset[f'PC{pc_y}'], 
                        label=f'Misclassified as {label}', 
                        alpha=0.6, edgecolor='red')

        plt.xlabel(f'PC{pc_x}')
        plt.ylabel(f'PC{pc_y}')
        plt.grid(True)

    plt.suptitle('Higher PC Comparisons: Correctly Classified vs Misclassified', fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Create a single legend for all subplots
    handles, labels = plt.gca().get_legend_handles_labels()
    unique_labels = dict(zip(labels, handles))
    plt.figlegend(unique_labels.values(), unique_labels.keys(), loc='upper center', ncol=3, fontsize='small', bbox_to_anchor=(0.5, 0.02))

    plt.savefig(f'{outpath}/compadre_1kg_pca_higher_pc_comparisons.png', bbox_inches='tight')
    plt.close()

    # Extract samples from all_pca that are not in kg_true IID column
    test_data = all_pca[~all_pca['IID'].isin(kg_true['IID'])]
    # Select only the numeric columns (e.g., PCs) for scaling
    test_data_numeric = test_data.iloc[:, 2:12].values
    test_data_scaled = scaler.transform(test_data_numeric)

    test_data['predicted_continental'] = clf.predict(test_data_scaled)
    y_proba = clf.predict_proba(test_data_scaled)

    # Plot the first two principal components for test data predictions
    plt.figure(figsize=(10, 8))

    # Plot each predicted population with a unique color
    predicted_populations = test_data['predicted_continental'].unique()
    for pop in predicted_populations:
        subset = test_data[test_data['predicted_continental'] == pop]
        plt.scatter(subset['PC1'], subset['PC2'], label=pop, alpha=0.6)

    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('Test Data PCA')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'{outpath}/compadre_test_data_pca_predictions.png')
    plt.close()

    # Plot higher PC comparisons for test data
    plt.figure(figsize=(12, 10))

    for i, (pc_x, pc_y) in enumerate(pairs, start=1):
        plt.subplot(2, 2, i)

        # Plot each predicted population with a unique color
        for pop in predicted_populations:
            subset = test_data[test_data['predicted_continental'] == pop]
            plt.scatter(subset[f'PC{pc_x}'], subset[f'PC{pc_y}'], label=pop, alpha=0.6)

        plt.xlabel(f'PC{pc_x}')
        plt.ylabel(f'PC{pc_y}')
        plt.grid(True)

    plt.suptitle('Higher PC Comparisons for Test Data', fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Create a single legend for all subplots
    handles, labels = plt.gca().get_legend_handles_labels()
    unique_labels = dict(zip(labels, handles))
    plt.figlegend(unique_labels.values(), unique_labels.keys(), loc='upper center', ncol=3, fontsize='small', bbox_to_anchor=(0.5, 0.02))

    plt.savefig(f'{outpath}/compadre_test_data_higher_pc_comparisons.png', bbox_inches='tight')
    plt.close()

    # Output predictions to a tab-separated table
    test_data = test_data.reset_index(drop=True)  # Reset index to ensure alignment
    output_file = f'{outpath}/compadre_test_data_predictions.tsv'
    with open(output_file, 'w') as f:
        f.write("SampleID\tPredictedPopulation\tClassificationProbabilities\n")
        for i, sample_id in enumerate(test_data['IID']):
            predicted_population = test_data.loc[i, 'predicted_continental']
            probabilities = ','.join(map(str, y_proba[i]))
            f.write(f"{sample_id}\t{predicted_population}\t{probabilities}\n")

    #print(f"Predictions saved to {output_file}")

#################################################################################

def run_new(pca_file, pop_file, verbose=False, linear=False):

    # Get prePRIMUS output directory to save figures 
    pca_file_bits = pca_file.split('/')[:-1]
    outpath = '/'.join(pca_file_bits)

    # Check if file has a header line (plink2 PCA output includes header, plink1.9 does not)
    # There's probably a more robust way to do this but for now it's fine

    with open(pca_file, 'r') as f:
        first_line = f.readline().strip()

    if first_line.startswith('#'): # Skip header
        all_pca = pd.read_csv(pca_file, header=None, sep='\s+', skiprows=1)
    else:
        all_pca = pd.read_csv(pca_file, header=None, sep='\s+')

    # Load PLINK PCA data
    num_pcs = all_pca.shape[1] - 2  # Subtract 2 for FID and IID columns
    #all_pca.columns = ['FID', 'IID'] + [f'PC{i}' for i in range(1, 21)]
    all_pca.columns = ['FID', 'IID'] + [f'PC{i}' for i in range(1, num_pcs + 1)]

    kg_true = pd.read_csv(pop_file, sep="\t", header=None)
    kg_true.columns = ['IID', 'POP']

    # Extract rows from all_pca where IID is in kg_true IID column
    kg_data = all_pca[all_pca['IID'].isin(kg_true['IID'])]
    kg_data = kg_data.merge(kg_true, on='IID')

    # Extract PCs 1-10
    pc_end_col = min(12, all_pca.shape[1])  # Use up to PC10 or whatever's available
    #X_train = kg_data.iloc[:, 2:12].values
    X_train = kg_data.iloc[:, 2:pc_end_col].values
    y_train = kg_data['POP']

    # Probably not necessary because of the way the data is prepared but better safe than sorry
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)

    # From testing these were the best performing settings using HM3 as the train set and 1kg as the test.
    if linear == False:
        clf = svm.SVC(kernel='rbf', C=5, gamma=0.5)
    else:
        clf = svm.LinearSVC(penalty='l2')
    clf = CalibratedClassifierCV(clf, cv=5)  # Use cross-validation for calibration
    clf.fit(X_train_scaled, y_train)

    # Perform 5-fold cross-validation and calculate accuracy
    cv_scores = cross_val_score(clf, X_train_scaled, y_train, cv=5, scoring='accuracy')
    if verbose == True:
        print(f"5-Fold Cross-Validation Accuracy: {np.mean(cv_scores):.4f} ± {np.std(cv_scores):.4f}")

    # Get true and predicted categories for kg_data
    kg_data['predicted'] = clf.predict(X_train_scaled)

    # Plot the first two principal components for 1000 Genomes data colored by population
    plt.figure(figsize=(10, 8))

    # Plot each population with a unique color
    populations = kg_data['POP'].unique()
    for pop in populations:
        subset = kg_data[kg_data['POP'] == pop]
        plt.scatter(subset['PC1'], subset['PC2'], label=pop, alpha=0.6)

    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('1000 Genomes PCA: Colored by Continental Population')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'{outpath}/compadre_PC_1kg_pca_by_population.png')
    plt.close()

    # Plot the first two principal components for 1000 Genomes data
    plt.figure(figsize=(10, 8))

    # Plot all correctly classified samples in grey
    correct_samples = kg_data[kg_data['POP'] == kg_data['predicted']]
    plt.scatter(correct_samples['PC1'], correct_samples['PC2'], label='Correctly Classified', color='grey', alpha=0.6)

    # Highlight misclassified samples with distinct colors
    misclassified_samples = kg_data[kg_data['POP'] != kg_data['predicted']]
    for label in np.unique(misclassified_samples['predicted']):
        subset = misclassified_samples[misclassified_samples['predicted'] == label]
        plt.scatter(subset['PC1'], subset['PC2'], label=f'Misclassified as {label}', alpha=0.6, edgecolor='red')

    # Annotate misclassified samples with the true population
    for _, row in misclassified_samples.iterrows():
        plt.annotate(row['POP'], (row['PC1'], row['PC2']), fontsize=8, color='black')

    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('1000 Genomes PCA: True Labels and Misclassifications')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True)
    plt.savefig(f'{outpath}/compadre_PC_1kg_pca_misclassifications.png', bbox_inches='tight')
    plt.close()

    # Plot higher PC comparisons (3:4, 5:6, etc.)
    plt.figure(figsize=(12, 10))

    pairs = [(3, 4), (5, 6), (7, 8), (9, 10)]

    for i, (pc_x, pc_y) in enumerate(pairs, start=1):
        plt.subplot(2, 2, i)

        # Plot each population with a unique color
        for pop in populations:
            subset = kg_data[kg_data['POP'] == pop]
            plt.scatter(subset[f'PC{pc_x}'], subset[f'PC{pc_y}'], label=pop, alpha=0.6)

        plt.xlabel(f'PC{pc_x}')
        plt.ylabel(f'PC{pc_y}')
        plt.grid(True)

    plt.suptitle('Higher PC Comparisons', fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Create a single legend for all subplots
    handles, labels = plt.gca().get_legend_handles_labels()
    unique_labels = dict(zip(labels, handles))
    plt.figlegend(unique_labels.values(), unique_labels.keys(), loc='upper center', ncol=3, fontsize='small', bbox_to_anchor=(0.5, 0.02))

    plt.savefig(f'{outpath}/compadre_PC_1kg_pca_higher_pc.png', bbox_inches='tight')
    plt.close()

    # Plot higher PC comparisons (3:4, 5:6, etc.) with misclassification
    plt.figure(figsize=(12, 10))

    for i, (pc_x, pc_y) in enumerate(pairs, start=1):
        plt.subplot(2, 2, i)

        # Plot all correctly classified samples in grey
        plt.scatter(correct_samples[f'PC{pc_x}'], correct_samples[f'PC{pc_y}'], 
                    label='Correctly Classified', 
                    color='grey', alpha=0.6)

        # Highlight misclassified samples with distinct colors
        for label in np.unique(misclassified_samples['predicted']):
            subset = misclassified_samples[misclassified_samples['predicted'] == label]
            plt.scatter(subset[f'PC{pc_x}'], subset[f'PC{pc_y}'], 
                        label=f'Misclassified as {label}', 
                        alpha=0.6, edgecolor='red')

        plt.xlabel(f'PC{pc_x}')
        plt.ylabel(f'PC{pc_y}')
        plt.grid(True)

    plt.suptitle('Higher PC Comparisons: Correctly Classified vs Misclassified', fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Create a single legend for all subplots
    handles, labels = plt.gca().get_legend_handles_labels()
    unique_labels = dict(zip(labels, handles))
    plt.figlegend(unique_labels.values(), unique_labels.keys(), loc='upper center', ncol=3, fontsize='small', bbox_to_anchor=(0.5, 0.02))

    plt.savefig(f'{outpath}/compadre_PC_1kg_pca_higher_pc_comparisons.png', bbox_inches='tight')
    plt.close()

    # Extract samples from all_pca that are not in kg_true IID column
    test_data = all_pca[~all_pca['IID'].isin(kg_true['IID'])]

    # Select only the numeric columns (e.g., PCs) for scaling
    #test_data_numeric = test_data.iloc[:, 2:12].values
    test_data_numeric = test_data.iloc[:, 2:pc_end_col].values
    test_data_scaled = scaler.transform(test_data_numeric)

    test_data['predicted'] = clf.predict(test_data_scaled)
    y_proba = clf.predict_proba(test_data_scaled)

    # Plot the first two principal components for test data predictions
    plt.figure(figsize=(10, 8))

    # Plot each predicted population with a unique color
    predicted_populations = test_data['predicted'].unique()
    for pop in predicted_populations:
        subset = test_data[test_data['predicted'] == pop]
        plt.scatter(subset['PC1'], subset['PC2'], label=pop, alpha=0.6)

    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('Test Data PCA')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'{outpath}/compadre_PC_pca_predictions.png')
    plt.close()

    # Plot higher PC comparisons for test data
    plt.figure(figsize=(12, 10))

    for i, (pc_x, pc_y) in enumerate(pairs, start=1):
        plt.subplot(2, 2, i)

        # Plot each predicted population with a unique color
        for pop in predicted_populations:
            subset = test_data[test_data['predicted'] == pop]
            plt.scatter(subset[f'PC{pc_x}'], subset[f'PC{pc_y}'], label=pop, alpha=0.6)

        plt.xlabel(f'PC{pc_x}')
        plt.ylabel(f'PC{pc_y}')
        plt.grid(True)

    plt.suptitle('Higher PC Comparisons', fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Create a single legend for all subplots
    handles, labels = plt.gca().get_legend_handles_labels()
    unique_labels = dict(zip(labels, handles))
    plt.figlegend(unique_labels.values(), unique_labels.keys(), loc='upper center', ncol=3, fontsize='small', bbox_to_anchor=(0.5, 0.02))

    plt.savefig(f'{outpath}/compadre_PC_higher_pc_comparisons.png', bbox_inches='tight')
    plt.close()

    # Output predictions to a tab-separated table
    test_data = test_data.reset_index(drop=True)  # Reset index to ensure alignment
    all_pops = []
    output_file = f'{outpath}/compadre_PC_population_predictions.tsv'
    with open(output_file, 'w') as f:
        f.write("SampleID\tPredictedPopulation\tClassificationProbabilities\n")
        for i, sample_id in enumerate(test_data['IID']):
            predicted_population = test_data.loc[i, 'predicted']
            all_pops.append(predicted_population)
            probabilities = ','.join(map(str, y_proba[i]))
            f.write(f"{sample_id}\t{predicted_population}\t{probabilities}\n")

    #print(f"Predictions saved to {output_file}")

    # RETURN PREDICTIONS HERE

    all_pops = list(set(all_pops))

    if verbose == True:
        print(f'ref_pops: {all_pops}\n')

    return all_pops


if __name__ == '__main__':

    #pca_file, pop_file = sys.argv[1], sys.argv[2]

    #pca_file_hm3 = '/data100t1/home/grahame/projects/compadre/unified-simulations/analysis-output/simulations_v3/EUR/uniform3_size20_sim100/uniform3_size20_sim100/primus/uniform3_size20_sim100_prePRIMUS/uniform3_size20_sim100_noDups_autosomal_merged.eigenvec'

    #pca_file_1kg = '/data100t1/home/grahame/projects/compadre/unified-simulations/analysis-output/simulations_v3/EUR/uniform3_size20_sim100/uniform3_size20_sim100/1kg_test/uniform3_size20_sim100_prePRIMUS/uniform3_size20_sim100_noDups_autosomal_merged.eigenvec'

    #pop_file = '/data100t1/home/grahame/projects/compadre/primus-ersa-v2/lib/1KG/1kg_phase3.id.txt'

    ############################

    #test1 = '/data100t1/home/grahame/projects/compadre/compadre-test/testing/full3/uniform3_size20_sim100_prePRIMUS/uniform3_size20_sim100_noDups_autosomal_merged.eigenvec'

    # AMR before and after LD pruning

    #test1 = '/data100t1/home/grahame/projects/compadre/compadre-test/testing/amr_full1/uniform3_size20_sim100_prePRIMUS/uniform3_size20_sim100_noDups_autosomal_merged.eigenvec'

    test1 = '/data100t1/home/grahame/projects/compadre/compadre-test/testing/amr_full3/uniform3_size20_sim100_prePRIMUS/uniform3_size20_sim100_noDups_autosomal_merged.eigenvec'

    # LD update 50 0.5 1
    test1 = '/data100t1/home/grahame/projects/compadre/compadre-test/testing/full1_newfolder2_LD/uniform3_size20_sim100_prePRIMUS/uniform3_size20_sim100_noDups_autosomal_merged.eigenvec'

    test2 = '/belowshare/vumcshare/data100t1/home/grahame/projects/compadre/compadre-test/lib/1KG/1KG_pop_classifier_ids.txt'

    run_new(test1, test2, verbose=True, linear=False)
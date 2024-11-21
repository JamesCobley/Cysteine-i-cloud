import pandas as pd

# List of Excel files to process
excel_files = [
    '/content/PTP1B_proteoform_final_distribution_with_counts (12).xlsx',
    '/content/PTP1B_proteoform_final_distribution_with_counts_1000nM (1).xlsx',
    '/content/PTP1B_proteoform_final_distribution_with_counts_100nM.xlsx'
]

# Load the CSV file containing the proteoform matrix
file_path_csv = '/content/PTP1B_proteoforms_ordered.csv'
df_matrix = pd.read_csv(file_path_csv)
print("Proteoform matrix loaded successfully.")

# Initialize a dictionary to store results for each task
results = {f"Task{i}": [] for i in range(1, 10)}

# Process each Excel file
for file_path_excel in excel_files:
    # Load the Excel file
    df_excel = pd.read_excel(file_path_excel, sheet_name='Sheet1')  # Adjust 'Sheet1' if needed
    print(f"Excel data from {file_path_excel} loaded successfully.")

    # Task 1: Define the number of proteoforms present and express as a percentage of 1,024
    num_proteoforms = df_excel['Proteoform_ID'].nunique()
    total_possible_proteoforms = 1024
    proteoform_percentage = (num_proteoforms / total_possible_proteoforms) * 100
    results['Task1'].append([file_path_excel, num_proteoforms, f"{proteoform_percentage:.2f}%"])

    # Task 2: Print the number of molecules in each k_value (write to table format)
    molecules_per_k_value = df_excel.groupby('k_value')['Count'].sum().reset_index()
    molecules_per_k_value['File'] = file_path_excel  # Add file origin as a column
    results['Task2'].append(molecules_per_k_value)

    # Task 3: Calculate the weighted mean redox state of the population based on k values (as a percentage)
    total_molecules = df_excel['Count'].sum()
    weighted_mean_redox = (df_excel['k_value'] * df_excel['Count']).sum() / total_molecules
    weighted_mean_redox_percentage = weighted_mean_redox * 10  # Convert to percentage
    results['Task3'].append([file_path_excel, f"{weighted_mean_redox_percentage:.2f}%"])

    # Task 4: Print the total number of molecules
    results['Task4'].append([file_path_excel, total_molecules])

    # Task 5: Print the number of molecules in proteoform "PF005" and calculate its percentage of the total
    proteoform_pf005 = df_excel[df_excel['Proteoform_ID'] == 'PF005']
    if not proteoform_pf005.empty:
        molecules_pf005 = proteoform_pf005['Count'].sum()
        pf005_percentage = (molecules_pf005 / total_molecules) * 100
        results['Task5'].append([file_path_excel, molecules_pf005, f"{pf005_percentage:.2f}%"])
    else:
        results['Task5'].append([file_path_excel, 'Proteoform PF005 not found', None])

    # Task 6: Calculate cysteine oxidation state percentages based on the enumerated matrix (write to table format)
    df_combined = pd.merge(df_excel, df_matrix, left_on='Proteoform_ID', right_on='Unnamed: 0')
    cysteine_columns = df_matrix.columns[1:]  # All cysteine columns (e.g., 'Cys32', 'Cys92', ...)
    oxidation_percentages = (df_combined[cysteine_columns].T * df_combined['Count']).sum(axis=1) / df_combined['Count'].sum() * 100

    # Create a DataFrame to store the oxidation percentages with the file name as a column
    oxidation_df = pd.DataFrame(oxidation_percentages).T
    oxidation_df.insert(0, 'File', file_path_excel)  # Insert file name as the first column
    results['Task6'].append(oxidation_df)

    # Task 7: Count proteoforms where Cys215 is oxidized (1)
    proteoforms_with_cys215_oxidized = df_combined[df_combined['Cys215'] == 1]
    num_proteoforms_cys215_oxidized = proteoforms_with_cys215_oxidized['Proteoform_ID'].nunique()
    percentage_cys215_oxidized = (num_proteoforms_cys215_oxidized / 1024) * 100
    results['Task7'].append([file_path_excel, num_proteoforms_cys215_oxidized, f"{percentage_cys215_oxidized:.2f}%"])

    # Task 8: Calculate percentage of PTP1B activation (proteoforms where Cys215 is oxidized as percentage of total population)
    molecules_with_cys215_oxidized = proteoforms_with_cys215_oxidized['Count'].sum()
    percentage_ptp1b_activation = (molecules_with_cys215_oxidized / total_molecules) * 100
    results['Task8'].append([file_path_excel, f"{percentage_ptp1b_activation:.2f}%"])

    # Task 9: Print the top 10 most frequent proteoform IDs and their k_values (write to table format)
    top_10_proteoforms = df_excel.nlargest(10, 'Count')[['Proteoform_ID', 'k_value']]
    top_10_proteoforms['File'] = file_path_excel  # Add file origin as a column
    results['Task9'].append(top_10_proteoforms)

# Write results to an Excel file with each task on a separate sheet
with pd.ExcelWriter("PTP1B_analysis_results.xlsx", engine='xlsxwriter') as writer:
    for task, data in results.items():
        # Handle DataFrame-type results (e.g., Task 2, Task 6, Task 9)
        if isinstance(data[0], pd.DataFrame):
            task_df = pd.concat(data, ignore_index=True)
            task_df.to_excel(writer, sheet_name=task, index=False)
        else:
            # Ensure all rows have 4 columns, adding None if necessary
            task_results = pd.DataFrame(data)
            task_results.to_excel(writer, sheet_name=task, index=False)

print("Analysis completed. Results saved to 'PTP1B_analysis_results.xlsx'.")

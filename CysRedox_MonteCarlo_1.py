import numpy as np
import pandas as pd
from google.cloud import storage
import time
from io import StringIO, BytesIO  # Added BytesIO for Excel file handling

# Google Cloud Storage setup
bucket_name = "jamesmontecarlo"
file_name = "PTP1B_proteoforms_ordered.csv"
storage_client = storage.Client()
bucket = storage_client.get_bucket(bucket_name)

# Load the proteoform library CSV file from Google Cloud Storage
blob = bucket.blob(file_name)
csv_data = blob.download_as_text()
proteoform_df = pd.read_csv(StringIO(csv_data))
proteoform_df['k_value'] = proteoform_df.iloc[:, 1:].sum(axis=1)

# Simulated parameters
number_of_molecules = 70000
time_steps = 300

# Provided oxidation and reduction probabilities
P_oxidation_Cys215 = 0.1254
P_oxidation_other = 0.0312
P_reduction_Cys215 = 0.3001
P_reduction_other = 0.1287

# Initial proteoform state for all molecules
initial_proteoform = 'PF001'

# Initialize a route map DataFrame
route_map = pd.DataFrame(index=np.arange(number_of_molecules), columns=np.arange(time_steps + 1))
route_map.iloc[:, 0] = initial_proteoform

# Initialize variables to track unique proteoforms and k distribution
unique_proteoforms = set()
final_k_distribution = np.zeros(11)
final_proteoform_counts = {}

# Utility function to find the proteoform ID based on state vector
def find_proteoform_id(state_vector):
    return proteoform_df.loc[(proteoform_df.iloc[:, 1:-1] == state_vector).all(axis=1), 'Unnamed: 0'].values[0]

# Utility function to check valid transitions
def valid_transition(current_state, target_state):
    diff = np.array(target_state) - np.array(current_state)
    oxidations = np.where(diff == 1)[0]
    reductions = np.where(diff == -1)[0]
    return len(oxidations) == 1 and len(reductions) == 0 or len(oxidations) == 0 and len(reductions) == 1

# Prepare proteoform dictionary for quick lookup
proteoform_library = proteoform_df.groupby(proteoform_df.iloc[:, 1:-1].sum(axis=1)).apply(lambda x: x.iloc[:, 1:-1].values.tolist()).to_dict()

# Progress function
def print_progress(step, total_steps):
    print(f"Step {step}/{total_steps} completed.")

# Simulation loop with time tracking
start_time = time.time()
for step in range(1, time_steps + 1):
    for i in range(number_of_molecules):
        current_state = route_map.iloc[i, step - 1]  
        current_vector = proteoform_df.loc[proteoform_df['Unnamed: 0'] == current_state].iloc[0, 1:-1].tolist()
        
        k_value = sum(current_vector)
        oxidized = False
        reduced = False
        new_vector = current_vector.copy()  

        # Oxidation step
        if k_value < max(proteoform_library.keys()):
            oxidation_prob = P_oxidation_Cys215 if current_vector[3] == 0 else P_oxidation_other
            if np.random.rand() < oxidation_prob:
                possible_states = proteoform_library[k_value + 1]
                valid_states = [s for s in possible_states if valid_transition(current_vector, s)]
                if valid_states:
                    new_vector = valid_states[np.random.randint(len(valid_states))]
                    new_proteoform = find_proteoform_id(new_vector)
                    oxidized = True

        # Reduction step
        if k_value > 0 and step > 1 and not oxidized:
            reduction_prob = P_reduction_Cys215 if current_vector[3] == 1 else P_reduction_other
            if np.random.rand() < reduction_prob:
                possible_states = proteoform_library[k_value - 1]
                valid_states = [s for s in possible_states if valid_transition(current_vector, s)]
                if valid_states:
                    new_vector = valid_states[np.random.randint(len(valid_states))]
                    new_proteoform = find_proteoform_id(new_vector)
                    reduced = True

        if not oxidized and not reduced:
            new_proteoform = current_state

        route_map.iloc[i, step] = new_proteoform

        unique_proteoforms.add(new_proteoform)
        
        if step == time_steps:
            final_k_distribution[sum(new_vector)] += 1
            if new_proteoform in final_proteoform_counts:
                final_proteoform_counts[new_proteoform] += 1
            else:
                final_proteoform_counts[new_proteoform] = 1

    print_progress(step, time_steps)

end_time = time.time()
print(f"Simulation completed in {end_time - start_time:.2f} seconds.")

# Calculate the overall redox state of the population
overall_redox_state = sum(k * count for k, count in enumerate(final_k_distribution)) / number_of_molecules

# Print the results
print(f"Total unique proteoforms formed: {len(unique_proteoforms)} out of 1,024")
print(f"Final k distribution:\n{final_k_distribution}")
print(f"Overall redox state of the population: {overall_redox_state}")

# Save the route map to a CSV file and upload to Google Cloud Storage
csv_buffer = StringIO()
route_map.to_csv(csv_buffer, index=False)
blob = bucket.blob("proteoform_route_map_POX_125_chaos.csv")
blob.upload_from_string(csv_buffer.getvalue())
print("Route map saved to Google Cloud Storage bucket 'jamesmontecarlo'")

# Create a DataFrame from the final proteoform counts and upload to Google Cloud Storage as Excel
final_counts_df = pd.DataFrame(list(final_proteoform_counts.items()), columns=['Proteoform_ID', 'Count'])
final_counts_df['k_value'] = final_counts_df['Proteoform_ID'].apply(lambda pid: sum(proteoform_df.loc[proteoform_df['Unnamed: 0'] == pid].iloc[0, 1:-1]))

# Save the results to Excel and upload to Google Cloud Storage
excel_buffer = BytesIO()  # Use BytesIO instead of StringIO
final_counts_df.to_excel(excel_buffer, index=False)

# Upload the Excel file to Google Cloud Storage
blob = bucket.blob("PTP1B_proteoform_final_distribution_with_counts_POX_125_chaos.xlsx")
blob.upload_from_string(excel_buffer.getvalue(), content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
print("Final distribution saved to Google Cloud Storage bucket 'jamesmontecarlo'")

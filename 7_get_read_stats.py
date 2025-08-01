import os
import glob
import subprocess
import pandas as pd


def count_reads_in_fastq(file_path):
    """Counts the number of reads in a FASTQ file using grep and wc."""
    result = subprocess.run(['grep', '-c', '^@', file_path], stdout=subprocess.PIPE, text=True)
    return int(result.stdout.strip())


def categorize_sample(file_name):
    """Categorizes the sample based on the file name."""
    if 'cow' in file_name:
        return 'cow'
    elif 'human' in file_name:
        return 'human'
    elif 'deer' in file_name:
        return 'deer'
    elif 'sheep' in file_name:
        return 'sheep'
    elif 'trim_fastp' in file_name:
        return 'virus'
    else:
        return 'unknown'


def main():
    # Define the directory containing the FASTQ files
    directory = '.'  # Update this to the correct path if needed

    # Initialize a list to hold the results
    results = []

    # Iterate through all FASTQ files in the directory
    for file_path in glob.glob(os.path.join(directory, "*.fastq")):
        file_name = os.path.basename(file_path)

        # Split the file name into parts
        parts = file_name.split('_')

        # Extract sample stem and read pair
        if len(parts) >= 2:
            sample_stem = '_'.join(parts[:-2])
            read_pair = parts[-2].replace('__1', 'R1').replace('__2', 'R2')  # Adjust for read pair naming
        else:
            sample_stem = parts[0]
            read_pair = "R1"  # Default to R1 if no explicit read pair is found

        # Count reads in the file
        reads_count = count_reads_in_fastq(file_path)

        # Categorize the sample
        category = categorize_sample(file_name)

        # Append the result
        results.append({
            'File Name': file_name,
            'Sample': sample_stem,
            'Category': category,
            'Read Pair': read_pair,
            'Read Count': reads_count
        })

    # Convert the results to a DataFrame
    df = pd.DataFrame(results)

    # Save the DataFrame to a CSV file
    df.to_csv('detailed_reads_counts.csv', index=False)

    # Print the resulting DataFrame
    print(df)


if __name__ == "__main__":
    main()

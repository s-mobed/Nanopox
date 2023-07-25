# Imports
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# Define a dictionary to convert units to bytes
byte_dict = {'B': 1, 'kB': 1000, 'MB': 1000000, 'GB': 1000000000, 'TB': 1000000000000}


def convert_to_bytes(value):
    if value != '0':
        number, unit = value.split()
        return float(number) * byte_dict[unit]
    else:
        return 0


def boxplot_scatter(dataframe, group_column, data_columns,
                    plot_x_label='', plot_y_label='', plot_title='', y_scale_down=1):

    # applies scale down
    dataframe[data_columns] = dataframe[data_columns] / y_scale_down

    # Group-by for boxplot
    grouped_data = dataframe.groupby(group_column)[data_columns]

    # saves boxplot data
    box_plot = plt.boxplot([grouped_data.get_group(process) for process in grouped_data.groups],
                           labels=grouped_data.groups.keys(), widths=0.3)

    # Plot elements
    plt.xlabel(plot_x_label)
    plt.ylabel(plot_y_label)
    plt.title(plot_title)
    plt.tight_layout()

    # Get the outliers for each process group, so to exclude for scatterplot
    outliers = [item.get_ydata() for item in box_plot['fliers']]

    # Plot the scatter points for non-outlier data
    for i, process in enumerate(grouped_data.groups):
        y_values = grouped_data.get_group(process)
        non_outlier_y = [y for y in y_values if y not in outliers[i]]
        x_values = np.random.normal(i + 1.3, 0.03, size=len(non_outlier_y))
        plt.scatter(x_values, non_outlier_y, edgecolor='black', facecolor='white', marker='o', alpha=0.8)

    # Display the plot
    plt.savefig(f"report_plots/{plot_title.replace(' ', '_')}.png")
    plt.close()


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Generate boxplot and scatter plots from a CSV file.")
    parser.add_argument("input_file", help="Path to the input CSV file")
    args = parser.parse_args()

    # Load CSV data and drop empty cells
    metrics = pd.read_csv(args.input_file).dropna()

    # Tidying the frame
    metrics['process'] = metrics['process'].str.split(':').str[1]

    # Processing the "duration" column to convert the seconds
    metrics['duration'] = metrics['duration'].apply(lambda x: '00:' + x if 'm' not in x else x).apply(
        lambda x: '00:' + x if 'h' not in x else x)
    metrics['duration'] = metrics['duration'].astype('str').str.replace(' ', '').str.replace('h', ':').str.replace('m',
                                                                                                                   ':').str.replace(
        's', '')
    metrics['duration'] = metrics['duration'].apply(
        lambda x: ':'.join(['0' + chunk if len(chunk) == 1 else chunk for chunk in x.split(':')]))
    metrics['duration'] = pd.to_timedelta(metrics['duration']).dt.total_seconds()

    # Standardising columns with different measures of data size

    # Extract the metric from the "vmem" column
    metric_series = metrics['vmem'].astype(str).str.split(' ').str[1]

    # Find the most common metric
    most_common_metric = metric_series.value_counts().idxmax()

    # Iterate over columns in the first row
    for column in metrics.columns:
        # Check if 'GB', 'MB', or 'kB' are present in the rows of the column and if the column is of object type
        if any(unit in str(value) for unit in ['GB', 'MB', 'kB'] for value in metrics[column]) and metrics[column].dtype == 'object':
            # Apply the conversion function to the column
            metrics[column] = metrics[column].apply(convert_to_bytes)


    # Create the "report_plots" directory if it doesn't exist
    os.makedirs("report_plots", exist_ok=True)

    # Generate plots
    boxplot_scatter(metrics, 'process', 'vmem', 'Processes', 'Memory Usage ('+most_common_metric+')', 'Total memory usage', byte_dict[most_common_metric])
    boxplot_scatter(metrics, 'process', '%cpu', 'Processes', '% single core CPU usage', 'CPU core usage by process')
    boxplot_scatter(metrics, 'process', 'duration', 'Processes', 'Process duration (min)', 'Process duration', 60)


if __name__ == "__main__":
    main()

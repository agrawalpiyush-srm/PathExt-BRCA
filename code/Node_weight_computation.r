# Modify the function to accept arguments from the command line
args <- commandArgs(trailingOnly = TRUE)

# The first argument will be the input file, and the second will be the output file
input_file <- args[1]
output_file <- args[2]

process_data <- function(input_file, output_file) {
  
  # Read input CSV file
  data <- read.csv(input_file, header = TRUE)
  
  # Extract columns
  gene <- data$gene
  control <- data$control
  case <- data$case
  
  # Calculate fold change and log fold change
  fc = (case / control)
  actual_fc = (case / control)
  fc[fc < 1] <- 1 / fc[fc < 1]  # Adjust for fold changes less than 1
  log_fc = log2(fc)
  log_fc = round(log_fc, digits = 5)
  
  # Create a new data frame with the original data, actual fold change, and log fold change
  new = data.frame(data, actual_fc, log_fc)
  
  # Convert to matrix and order based on a column (column 2 in this case)
  new_file = as.matrix(new)
  ordered = new_file[order(new_file[, 2]), ]
  
  # Convert back to data frame and apply smoothing
  aa = data.frame(ordered)
  fitting <- lowess(control, log_fc)
  bb = data.frame(aa$gene, fitting)
  
  # Create final data frame with additional columns
  new_datafile = data.frame(ordered, bb)
  colnames(new_datafile)[6] <- "sorted_gene"
  colnames(new_datafile)[7] <- "sorted_control"
  colnames(new_datafile)[8] <- "expected_logfc"
  
  # Calculate antilog and node weight
  antilog = abs(2 ^ new_datafile$expected_logfc)
  node_weight = as.numeric(new_datafile$actual_fc) / antilog
  
  # Create a final data frame with gene name and node_weight
  final_file = data.frame(gene = new_datafile$sorted_gene, node_weight = node_weight)
  
  # Write the final file to CSV
  write.csv(final_file, file = output_file, row.names = FALSE)
}

# Run the process_data function
process_data(input_file, output_file)

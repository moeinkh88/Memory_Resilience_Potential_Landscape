using XLSX
using CSV
using DataFrames

# File paths
xlsx_file = "Dataset_S1.xlsx"
csv_file = "sheet_Y.csv"

# Read sheet "Y" into a DataFrame
data = DataFrame(XLSX.readtable(xlsx_file, "Y"))

# Write to CSV
CSV.write(csv_file, data)

println("Saved sheet Y as sheet_Y.csv")


# Read sheet "MmuE"

# Open the workbook
xf = XLSX.readxlsx(xlsx_file)

# Access the raw sheet
sheet = xf["MmuE"]

# Convert to a 2D array of values
data = sheet[:]

# Suppose the interaction matrix is in rows 2:11, cols 2:11
interaction = DataFrame(data[3:13, 1:12], :auto)

# growth as a DataFrame
growth = DataFrame(Growth = data[3:13, 13])

# susceptibilities as a DataFrame
suscep = DataFrame(suscep = data[3:13, 14])


library(neighborweights)

# Check what functions are exported
print("Exported functions from neighborweights:")
print(ls("package:neighborweights"))

# Check if commute_time exists but is not exported
if (exists("commute_time", where = asNamespace("neighborweights"))) {
  print("\ncommute_time exists in neighborweights namespace but is not exported")
} else {
  print("\ncommute_time does not exist in neighborweights")
}

# Look for similar functions
exports <- ls("package:neighborweights")
commute_related <- grep("commute", exports, value = TRUE, ignore.case = TRUE)
if (length(commute_related) > 0) {
  print("\nFunctions related to 'commute':")
  print(commute_related)
}
library(gGnome)
library(data.table)
# Set the parent directory containing the subfolders

parent_dir <- "../../results/jabba/"
subfolders <- list.dirs(parent_dir, full.names = TRUE, recursive = FALSE)
for (subfolder in subfolders) {
  output_path <- file.path(subfolder, "events.csv")
  if (file.exists(output_path)) {
    message("Skipping ", subfolder, " — events.csv already exists.")
    next
  }
 jabba_path <- file.path(subfolder, "jabba.rds")
  if (file.exists(jabba_path)) {
    jabba <- readRDS(jabba_path)
    gg.jabba <- gG(jabba = jabba)
    events <- events(gg.jabba, verbose = TRUE)
    table <- as.data.frame(events$meta$events)
    fwrite(table, output_path)
    message("Processed and saved events for: ", subfolder)
  } else {
    message("No jabba.rds found in: ", subfolder)
  }
}

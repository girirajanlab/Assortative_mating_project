library(data.table)

# Pass filenames to function
filenames <- commandArgs(trailingOnly=TRUE)
if (length(filenames)!=3) {
	stop('You must provide 3 file names to this function - QC file, valid file, and invalid file')
}

QC_file<- filenames[1]
valid_file <- filenames[2]
invalid_file <- filenames[3]

dat <- fread(QC_file)

# Get samples with F coefficient within 2SD of the population mean
valid <- dat[F<=mean(F)+2*sd(F) & F>=mean(F)-2*sd(F)]
invalid <- dat[F>mean(F)+2*sd(F) | F<mean(F)-2*sd(F)]

# print FID and IID for valid samples
fwrite(valid[,c("FID","IID")], valid_file, sep="\t")
fwrite(invalid[,c("FID","IID","F")], invalid_file, sep="\t")

# print number of samples that passed/failed QC
cat(paste(length(rownames(dat)), ' samples input.\n', length(rownames(valid)), ' passed heterozygosity filter.\n', length(rownames(invalid)), ' failed heterozygosity filter.\n', sep=''))

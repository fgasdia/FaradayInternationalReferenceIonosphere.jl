#==
!!! note

    This code does not need to be run by normal users.

This script converts the FIRI-2018 xlsx from https://figshare.com/s/357cb03b3e5bed649bbc
to CSV format.

This new file is then put into a .tar.gz and uploaded to my figshare account to be used
as a source URL for artifact. The original data is under a CC BY 4.0 license which gives
freedom to share (copy and redistribute in any medium or format) and adapt as long as
appropriate credit is provided.

Finally, ArtifactUtils.jl is used to generate the Artifact.toml file for
FaradayInternationalReferenceIonosphere.jl.
==#

# These packages are not all in FaradayInternationalReferenceIonosphere Project.toml
using Artifacts, ArtifactUtils
using CSV, DataFrames, XLSX
using Tar, CodecZlib

# User friendly URL is: https://figshare.com/s/357cb03b3e5bed649bbc
firi2018_url = "https://ndownloader.figshare.com/files/11823206?private_link=357cb03b3e5bed649bbc"

firi_xlsx = "firi2018.xlsx"
firi_targz = "firi.tar.gz"
firi_dir = "firi"
isdir(firi_dir) || mkdir(firi_dir)

download("$(firi2018_url)", firi_xlsx)

# Reading the XLSX is very slow
df = DataFrame(XLSX.readtable(firi_xlsx, "FIRI-new_vs_altitude")...)

CSV.write(joinpath(firi_dir, "firi2018.csv"), df)

# Tarball the csv
targz = open(firi_targz, "w")
tar = GzipCompressorStream(targz)
Tar.create(firi_dir, tar)
close(tar)

# Cleanup
rm(firi_xlsx)
rm(firi_dir, recursive=true)

# Now post firi.tar.gz to figshare.
# Friendly URL: https://figshare.com/s/33a146841a3f74a74590
# Then (don't skip posting to figshare, because it downloads the tarball!)...

# Write the Artifacts.toml file
add_artifact!(
    "Artifacts.toml",
    "firi",
    "https://ndownloader.figshare.com/files/26640602?private_link=33a146841a3f74a74590",
    force=true,
)

rm(firi_targz)
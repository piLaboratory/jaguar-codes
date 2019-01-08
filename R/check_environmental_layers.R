#-------------------------------------------------------------------
#
# Verifying spatial variables for the Neotropical Jaguar Movement Projects
#
# Bernardo Niebuhr
# Dec. 2018
#--------------------------------------------------------------------

# Load packages

# Root folder
wd <- 'H:/_neojaguardatabase/Envdatabase/Basesfinais_res30m'

setwd(wd)

# List individual folders
# Each folder has the environmental layers for the urroundings of one individual:
# 70 km buffer around the jaguar locations

dirs <- ?list.dirs('.')
(dirs <- dirs[-1]) # remove the root folder ('.')
length(dirs) # ok, 117 individuals

# List files per folder
li <- list()

for(i in 1:length(dirs)) {
  fil <- list.files(dirs[i], pattern = '.tif')
  li[[i]] <- fil[endsWith(fil, '.tif')]
}

# Check number of layers per folder
length(li)

# Automatized process for all individuals
sapply(li, length)
all(sapply(li, length) == length(li[[1]])) # ok

# Check layer names consistency

# Compare variables names, pairwise, ignoring the first characters that represent the individual ID
all(substring(li[[1]], 5) == substring(li[[2]], 5)) 

# Automatized process for all individuals
sapply(li, function(x) all(substring(x, 5) == substring(li[[1]], 5)))
all(sapply(li, function(x) all(substring(x, 5) == substring(li[[1]], 5)))) # ok!


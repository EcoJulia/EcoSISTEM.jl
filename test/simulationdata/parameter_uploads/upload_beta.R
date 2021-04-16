#'
#' Upload script for beta parameters - to control the impact of force of infection
#' and environmental reservoir
#'

library(SCRCdataAPI)

# Go to data.scrc.uk, click on Links, then Generate API Token, and save your
# token in your working directory as token.txt. If the following returns an
# error, then save a carriage return after the token.
key <- readLines("token.txt")
namespace <- "Simulation"

# Now I'm assuming you've already made a multi-parameter toml. Let's read it
# into the global envirnment so we can extract component_names from it
filename <- "beta.toml"
path <- "data-raw"
dat <- configr::read.config(file.path(path, filename))
component_names <- names(dat)

# The product_name is used to identify the data product and will be used to
# generate the file location in the ScottishCovidResponse/DataRepository
# GitHub repository:
product_name <- paste0("virus/beta")

# Remember to push the toml file to
# ScottishCovidResponse/DataRepository/[namespace]/[product_name]
# (e.g. ScottishCovidResponse/DataRepository/SCRC/human/infection/SARS-CoV-2)
# after completing this process

# The version number of the data product
productVersion <- "0.1.0"

# *******************************************************************
# Now run the code below (you probably don't need to change anything)
# *******************************************************************



# default data that should be in database ---------------------------------

# Assuming the toml will be stored in the ScottishCovidResponse/DataRepository
# GitHub repository
productStorageRoot <- "DataRepository"

storage_rootId <- new_storage_root(
  name = productStorageRoot,
  root = "https://raw.githubusercontent.com/ScottishCovidResponse/DataRepository/",
  key = key)

namespaceId <- new_namespace(name = namespace,
                             key = key)


# upload data product metadata to database --------------------------------

product_path <- paste("master", namespace, product_name, sep = "/")

upload_data_product(storage_root_id = storage_rootId,
                    name = product_name,
                    component_name = component_names,
                    processed_path = file.path(path, filename), # local
                    product_path = file.path(product_path, filename), # ftp
                    version = productVersion,
                    namespace_id = namespaceId,
                    key = key)


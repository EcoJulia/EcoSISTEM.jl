#' Virus decay rate
#'
#' Add a point-estimate to the data registry
#'

library(SCRCdataAPI)

# Go to data.scrc.uk, click on Links, then Generate API Token, and save your
# token in your working directory as token.txt. If the following returns an
# error, then save a carriage return after the token.
key <- readLines("token.txt")
namespace <- "Simulation"

# The product_name is used to identify the data product and will be used to
# generate the file location in the ScottishCovidResponse/DataRepository
# GitHub repository. In addition to this, when you create your toml file
# below, it will be saved to data-raw/[product_name].
product_name <- "virus/virus-decay"

# Remember to push the toml file to
# ScottishCovidResponse/DataRepository/[namespace]/[product_name]
# (e.g. ScottishCovidResponse/DataRepository/SCRC/human/infection/SARS-CoV-2)
# after completing this process

# The component_name is taken as the last part of the product_name
component_name <- gsub("^.*/([^/]*)$", "\\1", product_name)

# The value of the point-estimate
component_value <- 0.333
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


# generate toml -----------------------------------------------------------

path <- paste("master", namespace, product_name, sep = "/")
filename <- paste0(productVersion, ".toml")

create_estimate(filename = filename,
                path = file.path("data-raw", path),
                parameters = as.list(setNames(component_value, component_name)))


# upload data product metadata to database --------------------------------

upload_data_product(storage_root_id = storage_rootId,
                    name = product_name,
                    component_name = component_name,
                    processed_path = file.path("data-raw", path, filename),
                    product_path = file.path(path, filename),
                    version = productVersion,
                    namespace_id = namespaceId,
                    key = key)

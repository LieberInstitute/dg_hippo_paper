library('readr')
library('purrr')
library('dplyr')
library('stringr')
library('testthat')
library('sessioninfo')

## Main options to work through
regions <- c('DentateGyrus')
types <- c('psycm')
features <- c('gene', 'exon', 'jxn', 'tx')
file_types <- c('all', 'included', 'dropped')

## Function for locating the files
locate_files <- function(ftype, path, type) {
    
    ## Determine the file pattern to use
    fpatt <- case_when(
        ftype == 'all' ~ '',
        ftype == 'included' ~ '\\.analysis\\.joint_included',
        ftype == 'dropped' ~ '\\.analysis\\.joint_dropped'
    )
    patt <- paste0(
        type,
        '\\.[[:digit:]]*',
        fpatt, 
        '\\.dat$'
    )

    ## Find the files
    dat_files <- dir(
        path = path,
        pattern = patt,
        full.names = TRUE
    )
    
    ## Extract the chromosome
    names(dat_files) <- str_extract(
        basename(dat_files),
        '[:digit:]+'
    )
    
    ## Done
    return(dat_files)
}

## Check that the file locator is working
test_files <- locate_files(
    'all',
    '/dcl01/ajaffe/data/lab/dg_hippo_paper/twas/DentateGyrus/gene/psycm',
    'psycm'
)

test_that('File locator', {
    expect_equal(
        basename(test_files['11']),
        'psycm.11.dat'
    )
    expect_equivalent(
        test_files['1'],
        '/dcl01/ajaffe/data/lab/dg_hippo_paper/twas/DentateGyrus/gene/psycm/psycm.1.dat'
    )
})


## Now read in the data for all file types
## Note that each file type has different numbers of columns
## when why I'm keeping them as different elements of the
## 'twas' list
twas <- map(file_types, function(ftype) {
    
    ## data.frame with all the combinations of arguments
    arg_grid <- expand.grid(
        region = regions,
        type = types,
        feature = features,
        stringsAsFactors = FALSE
    )
    
    
    pmap_dfr(arg_grid, function(region, type, feature) {
        
        ## Construct the path to the files
        path <- file.path(
            '/dcl01/ajaffe/data/lab/dg_hippo_paper/twas',
            region,
            feature,
            type
        )
        
        ## Now locate the files
        dat_files <- locate_files(
            ftype = ftype,
            path = path,
            type = type
        )
        
        ## Next read the files and add the chromosome info
        result <- map2_dfr(
            dat_files,
            names(dat_files),
            function(f, chr) {
                res <- suppressMessages(read_tsv(f))
                res$chr <- chr
                return(res)
            }
        )
        
        ## Next add the region, feature and type information
        result$region <- region
        result$feature <- feature
        result$type <- type
        
        ## Done
        return(result)
    })
})
names(twas) <- file_types

## Explore the resulting dimensions
map_dfr(twas, dim)
# # A tibble: 2 x 3
#      all included dropped
#    <int>    <int>   <int>
# 1 68199      375   40654
# 2    24       12      12


## Save the data for later use
save(twas, file = 'rdas/twas.Rdata')

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

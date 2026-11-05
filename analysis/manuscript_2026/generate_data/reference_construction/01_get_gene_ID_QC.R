############################################################
# get_gene_ID.R
#
# Purpose:
#   Human gene symbol list
#       -> HGNC official symbol exact matching
#       -> NCBI Gene ID
#       -> NCBI ESummary reciprocal validation
#       -> final Gene ID list
#
# IMPORTANT:
#   Final Gene ID list is generated ONLY when all genes pass QC.
############################################################


rm(list = ls())

library(dplyr)
library(httr)
library(jsonlite)


############################################################
# 0. Settings
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop(
    paste0(
      "\nUsage:\n",
      "Rscript get_gene_ID.R ",
      "<ISG_symbol_list> ",
      "<Control_symbol_list> ",
      "<NCBI_email> ",
      "<output_geneID_list>\n"
    )
  )
}

CORE_FILE    <- args[1]
CONTROL_FILE <- args[2]
NCBI_EMAIL   <- args[3]
OUT_FINAL    <- args[4]

# Fixed reference
HGNC_FILE <- "HGNC_Ens_Ref_mod.list"

# Automatically generated QC files
out_prefix <- tools::file_path_sans_ext(OUT_FINAL)
OUT_QC      <- paste0(out_prefix, "_QC.tsv")
OUT_PROBLEM <- paste0(out_prefix, "_problem.tsv")

# Fixed settings
NCBI_TOOL <- "ISG_profiler_db_builder"
HUMAN_TAXID <- "9606"

############################################################
# 1. Function:
#    NCBI Gene ID -> NCBI Gene information
#
# Actual NCBI ESummary JSON structure confirmed:
#
# result
#   └── "103"
#        ├── uid
#        ├── name
#        ├── nomenclaturesymbol
#        ├── nomenclaturestatus
#        └── organism
#             └── taxid
############################################################

get_ncbi_gene_info <- function(
    gene_id,
    max_retry = 5
) {

  base_url <-
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


  for (attempt in seq_len(max_retry)) {

    result <- tryCatch({

      ########################################################
      # Query NCBI
      ########################################################

      res <- GET(
        base_url,
        query = list(
          db = "gene",
          id = as.character(gene_id),
          retmode = "json",
          tool = NCBI_TOOL,
          email = NCBI_EMAIL
        ),
        timeout(30)
      )


      stop_for_status(res)


      ########################################################
      # Parse JSON
      ########################################################

      txt <- content(
        res,
        as = "text",
        encoding = "UTF-8"
      )


      js <- fromJSON(
        txt,
        simplifyVector = FALSE
      )


      ########################################################
      # Get returned UID
      ########################################################

      if (
        is.null(js$result$uids) ||
        length(js$result$uids) == 0
      ) {

        stop(
          "No UID returned from NCBI."
        )
      }


      uid <- as.character(
        js$result$uids[[1]]
      )


      ########################################################
      # Get Gene record
      ########################################################

      rec <- js$result[[uid]]


      if (is.null(rec)) {

        stop(
          "Gene record was not found in NCBI response."
        )
      }


      ########################################################
      # Extract fields safely
      ########################################################

      ncbi_name <- NA_character_

      if (
        !is.null(rec$name) &&
        length(rec$name) > 0 &&
        rec$name != ""
      ) {

        ncbi_name <-
          as.character(rec$name)
      }


      ncbi_symbol <- NA_character_

      if (
        !is.null(rec$nomenclaturesymbol) &&
        length(rec$nomenclaturesymbol) > 0 &&
        rec$nomenclaturesymbol != ""
      ) {

        ncbi_symbol <-
          as.character(
            rec$nomenclaturesymbol
          )
      }


      ncbi_nomenclature_status <-
        NA_character_

      if (
        !is.null(rec$nomenclaturestatus) &&
        length(rec$nomenclaturestatus) > 0 &&
        rec$nomenclaturestatus != ""
      ) {

        ncbi_nomenclature_status <-
          as.character(
            rec$nomenclaturestatus
          )
      }


      ncbi_taxid <- NA_character_

      if (
        !is.null(rec$organism) &&
        !is.null(rec$organism$taxid)
      ) {

        ncbi_taxid <-
          as.character(
            rec$organism$taxid
          )
      }


      ########################################################
      # Return one-row data.frame
      ########################################################

      data.frame(

        NCBI_geneID =
          as.character(uid),

        NCBI_name =
          ncbi_name,

        NCBI_symbol =
          ncbi_symbol,

        NCBI_nomenclature_status =
          ncbi_nomenclature_status,

        NCBI_taxID =
          ncbi_taxid,

        NCBI_request_status =
          "OK",

        stringsAsFactors = FALSE
      )

    },


    ##########################################################
    # Error handling
    ##########################################################

    error = function(e) {

      message(
        "NCBI request failed for Gene ID ",
        gene_id,
        " (attempt ",
        attempt,
        "/",
        max_retry,
        "): ",
        conditionMessage(e)
      )

      NULL
    })


    ##########################################################
    # Successful request
    ##########################################################

    if (!is.null(result)) {

      return(result)
    }


    ##########################################################
    # Retry wait
    ##########################################################

    Sys.sleep(
      5 * attempt
    )
  }


  ############################################################
  # All attempts failed
  ############################################################

  data.frame(

    NCBI_geneID =
      as.character(gene_id),

    NCBI_name =
      NA_character_,

    NCBI_symbol =
      NA_character_,

    NCBI_nomenclature_status =
      NA_character_,

    NCBI_taxID =
      NA_character_,

    NCBI_request_status =
      "REQUEST_FAILED",

    stringsAsFactors = FALSE
  )
}



############################################################
# 2. Read HGNC reference
############################################################

cat("\n")
cat("========================================\n")
cat("Reading HGNC reference\n")
cat("========================================\n")


hgnc <- read.delim(
  HGNC_FILE,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


############################################################
# Required columns
############################################################

required_hgnc_columns <-
  c(
    "hum_symbol_HGNC",
    "hum_geneID"
  )


missing_columns <-
  setdiff(
    required_hgnc_columns,
    colnames(hgnc)
  )


if (length(missing_columns) > 0) {

  stop(
    paste(
      "Missing required HGNC columns:",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}


############################################################
# Normalize columns
############################################################

hgnc <- hgnc %>%
  mutate(

    hum_symbol_HGNC =
      trimws(
        as.character(
          hum_symbol_HGNC
        )
      ),

    hum_geneID =
      trimws(
        as.character(
          hum_geneID
        )
      )
  )


cat(
  "HGNC reference rows:",
  nrow(hgnc),
  "\n"
)



############################################################
# 3. HGNC reference QC
############################################################

cat("\n")
cat("========================================\n")
cat("HGNC reference QC\n")
cat("========================================\n")


############################################################
# Missing symbols
############################################################

missing_hgnc_symbol <-
  hgnc %>%
  filter(
    is.na(hum_symbol_HGNC) |
      hum_symbol_HGNC == ""
  )


cat(
  "Missing HGNC symbols:",
  nrow(missing_hgnc_symbol),
  "\n"
)


############################################################
# Duplicate official symbols
############################################################

dup_hgnc_symbol <-
  hgnc %>%
  filter(
    !is.na(hum_symbol_HGNC),
    hum_symbol_HGNC != ""
  ) %>%
  count(
    hum_symbol_HGNC,
    name = "n"
  ) %>%
  filter(
    n > 1
  )


cat(
  "Duplicated HGNC symbols:",
  nrow(dup_hgnc_symbol),
  "\n"
)


if (nrow(dup_hgnc_symbol) > 0) {

  print(
    as.data.frame(
      dup_hgnc_symbol
    )
  )

  stop(
    "Duplicated official symbols detected in HGNC reference."
  )
}



############################################################
# 4. Read Core ISG input
############################################################

core <- read.table(
  CORE_FILE,
  header = FALSE,
  stringsAsFactors = FALSE,
  fill = TRUE,
  comment.char = "",
  quote = ""
)


if (ncol(core) != 1) {

  stop(
    "CORE_FILE must contain exactly one column of gene symbols."
  )
}


colnames(core) <-
  "hum_symbol"


core <- core %>%
  transmute(

    hum_symbol =
      trimws(
        as.character(
          hum_symbol
        )
      ),

    type =
      "Core_ISG"
  ) %>%
  filter(
    !is.na(hum_symbol),
    hum_symbol != ""
  )



############################################################
# 5. Read Control input
############################################################

control <- read.table(
  CONTROL_FILE,
  header = FALSE,
  stringsAsFactors = FALSE,
  fill = TRUE,
  comment.char = "",
  quote = ""
)


if (ncol(control) != 1) {

  stop(
    "CONTROL_FILE must contain exactly one column of gene symbols."
  )
}


colnames(control) <-
  "hum_symbol"


control <- control %>%
  transmute(

    hum_symbol =
      trimws(
        as.character(
          hum_symbol
        )
      ),

    type =
      "Control"
  ) %>%
  filter(
    !is.na(hum_symbol),
    hum_symbol != ""
  )



############################################################
# 6. Input QC
############################################################

cat("\n")
cat("========================================\n")
cat("Input gene symbol QC\n")
cat("========================================\n")


cat(
  "Core ISGs:",
  nrow(core),
  "\n"
)


cat(
  "Controls:",
  nrow(control),
  "\n"
)

############################################################
# Duplicate Core symbols
############################################################

dup_core <-
  core %>%
  count(
    hum_symbol,
    name = "n"
  ) %>%
  filter(
    n > 1
  )


cat(
  "Duplicated Core symbols:",
  nrow(dup_core),
  "\n"
)


if (nrow(dup_core) > 0) {

  print(
    as.data.frame(
      dup_core
    )
  )
}



############################################################
# Duplicate Control symbols
############################################################

dup_control <-
  control %>%
  count(
    hum_symbol,
    name = "n"
  ) %>%
  filter(
    n > 1
  )


cat(
  "Duplicated Control symbols:",
  nrow(dup_control),
  "\n"
)


if (nrow(dup_control) > 0) {

  print(
    as.data.frame(
      dup_control
    )
  )
}



############################################################
# Core / Control overlap
############################################################

overlap <-
  intersect(
    core$hum_symbol,
    control$hum_symbol
  )


cat(
  "Core / Control overlap:",
  length(overlap),
  "\n"
)


if (length(overlap) > 0) {

  print(overlap)
}



############################################################
# Stop if input is structurally invalid
############################################################

if (
  nrow(dup_core) > 0 ||
  nrow(dup_control) > 0 ||
  length(overlap) > 0
) {

  stop(
    "Input gene symbol list failed QC."
  )
}



############################################################
# 7. Combine input genes
############################################################

genes <-
  bind_rows(
    core,
    control
  )


cat(
  "Total input genes:",
  nrow(genes),
  "\n"
)


cat(
  "Unique input symbols:",
  n_distinct(
    genes$hum_symbol
  ),
  "\n"
)



############################################################
# 8. HGNC official-symbol exact matching
############################################################

cat("\n")
cat("========================================\n")
cat("HGNC official-symbol mapping\n")
cat("========================================\n")


gene_map <-
  genes %>%
  left_join(

    hgnc %>%
      select(
        hum_symbol_HGNC,
        hum_geneID
      ),

    by = c(
      "hum_symbol" =
        "hum_symbol_HGNC"
    )
  )



############################################################
# 9. HGNC mapping QC
############################################################

gene_map <-
  gene_map %>%
  mutate(

    HGNC_QC =
      case_when(

        is.na(hum_geneID) |
          hum_geneID == "" ~
          "NOT_FOUND",

        TRUE ~
          "PASS"
      )
  )


cat("\nHGNC mapping summary:\n")


print(
  as.data.frame(
    gene_map %>%
      count(
        type,
        HGNC_QC
      )
  )
)



############################################################
# HGNC mapping problems
############################################################

hgnc_problem <-
  gene_map %>%
  filter(
    HGNC_QC != "PASS"
  )


if (nrow(hgnc_problem) > 0) {

  cat("\n")
  cat("=== HGNC mapping problems ===\n")

  print(
    as.data.frame(
      hgnc_problem
    )
  )


  write.table(
    hgnc_problem,
    OUT_PROBLEM,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )


  stop(
    "HGNC mapping failed. NCBI QC was not performed."
  )
}



############################################################
# 10. Duplicated Gene ID QC
############################################################

dup_geneID <-
  gene_map %>%
  count(
    hum_geneID,
    name = "n"
  ) %>%
  filter(
    n > 1
  )


cat(
  "\nDuplicated mapped Gene IDs:",
  nrow(dup_geneID),
  "\n"
)


if (nrow(dup_geneID) > 0) {

  print(
    as.data.frame(
      dup_geneID
    )
  )


  stop(
    "Multiple input symbols mapped to the same Gene ID."
  )
}



############################################################
# 11. NCBI reciprocal validation
############################################################

cat("\n")
cat("========================================\n")
cat("NCBI reciprocal Gene ID validation\n")
cat("========================================\n")


unique_geneIDs <-
  unique(
    gene_map$hum_geneID
  )


ncbi_results <-
  vector(
    "list",
    length(
      unique_geneIDs
    )
  )


for (
  i in seq_along(
    unique_geneIDs
  )
) {

  gene_id <-
    unique_geneIDs[i]


  cat(
    "[",
    i,
    "/",
    length(unique_geneIDs),
    "] Gene ID ",
    gene_id,
    "\n",
    sep = ""
  )


  ncbi_results[[i]] <-
    get_ncbi_gene_info(
      gene_id
    )


  # Keep requests conservative without an NCBI API key.
  Sys.sleep(0.4)
}


ncbi_results <-
  bind_rows(
    ncbi_results
  )



############################################################
# 12. Preserve returned NCBI Gene ID
############################################################

ncbi_results <-
  ncbi_results %>%
  mutate(

    NCBI_returned_geneID =
      as.character(
        NCBI_geneID
      )
  )



############################################################
# 13. Merge NCBI results with HGNC mapping
############################################################

gene_map <-
  gene_map %>%
  mutate(

    hum_geneID =
      as.character(
        hum_geneID
      )
  ) %>%

  left_join(

    ncbi_results,

    by = c(
      "hum_geneID" =
        "NCBI_geneID"
    )
  )



############################################################
# 14. Sanity check after NCBI merge
############################################################

cat("\n")
cat("========================================\n")
cat("NCBI returned data: first 10 genes\n")
cat("========================================\n")


display_columns <-
  c(
    "hum_symbol",
    "hum_geneID",
    "NCBI_returned_geneID",
    "NCBI_name",
    "NCBI_symbol",
    "NCBI_nomenclature_status",
    "NCBI_taxID",
    "NCBI_request_status"
  )


print(
  as.data.frame(
    head(
      gene_map[
        ,
        display_columns,
        drop = FALSE
      ],
      10
    )
  )
)



############################################################
# 15. NCBI reciprocal QC
#
# Check:
#
# 1. NCBI request succeeded
# 2. NCBI returned a Gene ID
# 3. HGNC Gene ID == NCBI returned Gene ID
# 4. Organism == Homo sapiens (taxid 9606)
# 5. NCBI has an official nomenclature symbol
# 6. Nomenclature status == Official
# 7. Input symbol == NCBI official symbol
############################################################

gene_map <-
  gene_map %>%
  mutate(

    NCBI_QC =
      case_when(

        is.na(NCBI_request_status) ~
          "NCBI_NOT_FOUND",

        NCBI_request_status ==
          "REQUEST_FAILED" ~
          "REQUEST_FAILED",

        is.na(NCBI_returned_geneID) |
          NCBI_returned_geneID == "" ~
          "NCBI_NOT_FOUND",

        as.character(hum_geneID) !=
          as.character(
            NCBI_returned_geneID
          ) ~
          "GENE_ID_MISMATCH",

        is.na(NCBI_taxID) |
          NCBI_taxID == "" ~
          "NO_TAXID",

        as.character(NCBI_taxID) !=
          HUMAN_TAXID ~
          "WRONG_TAXID",

        is.na(NCBI_symbol) |
          NCBI_symbol == "" ~
          "NO_OFFICIAL_SYMBOL",

        is.na(
          NCBI_nomenclature_status
        ) |
          NCBI_nomenclature_status == "" ~
          "NO_NOMENCLATURE_STATUS",

        NCBI_nomenclature_status !=
          "Official" ~
          "NOT_OFFICIAL",

        hum_symbol !=
          NCBI_symbol ~
          "SYMBOL_MISMATCH",

        TRUE ~
          "PASS"
      ),


    FINAL_QC =
      case_when(

        HGNC_QC == "PASS" &
          NCBI_QC == "PASS" ~
          "PASS",

        TRUE ~
          "FAIL"
      )
  )



############################################################
# 16. QC summary
############################################################

cat("\n")
cat("========================================\n")
cat("NCBI QC summary\n")
cat("========================================\n")


print(
  as.data.frame(
    gene_map %>%
      count(
        NCBI_QC
      )
  )
)


cat("\n")
cat("Final QC summary:\n")


print(
  as.data.frame(
    gene_map %>%
      count(
        type,
        FINAL_QC
      )
  )
)



############################################################
# 17. Identify problematic genes
############################################################

problem <-
  gene_map %>%
  filter(
    FINAL_QC != "PASS"
  )


cat("\n")
cat("=== Genes requiring review ===\n")


if (nrow(problem) == 0) {

  cat("None\n")

} else {

  print(
    as.data.frame(
      problem
    )
  )
}



############################################################
# 18. Save full QC results
############################################################

write.table(
  gene_map,
  OUT_QC,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  problem,
  OUT_PROBLEM,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat(
  "\nFull QC table saved:",
  OUT_QC,
  "\n"
)


cat(
  "Problem table saved:",
  OUT_PROBLEM,
  "\n"
)



############################################################
# 19. FAIL-FAST
#
# IMPORTANT:
# Final Gene ID list must NOT be generated if any gene fails.
############################################################

if (nrow(problem) > 0) {

  cat("\n")
  cat("========================================\n")
  cat("QC FAILED\n")
  cat("========================================\n")


  cat(
    "Problematic genes:",
    nrow(problem),
    "\n"
  )


  cat(
    "Final Gene ID list was NOT generated.\n"
  )


  stop(
    "HGNC/NCBI validation failed."
  )
}



############################################################
# 20. Create final Gene ID list
#
# REQUIRED OUTPUT FORMAT:
#
# hum_geneID    type
#
# No header
#
# Core_ISG -> ISG
############################################################

final_gene_list <-
  gene_map %>%
  transmute(

    hum_geneID =
      hum_geneID,

    type =
      if_else(
        type == "Core_ISG",
        "ISG",
        type
      )
  )



############################################################
# 21. Final output QC
############################################################

cat("\n")
cat("========================================\n")
cat("Final Gene ID list QC\n")
cat("========================================\n")


cat(
  "Rows:",
  nrow(final_gene_list),
  "\n"
)


cat(
  "Unique Gene IDs:",
  n_distinct(
    final_gene_list$hum_geneID
  ),
  "\n"
)


print(
  as.data.frame(
    final_gene_list %>%
      count(
        type
      )
  )
)



############################################################
# Defensive final checks
############################################################

if (
  nrow(final_gene_list) !=
    nrow(genes)
) {

  stop(
    "Final output row count does not match input gene count."
  )
}


if (
  n_distinct(
    final_gene_list$hum_geneID
  ) !=
    nrow(final_gene_list)
) {

  stop(
    "Final Gene ID list contains duplicated Gene IDs."
  )
}



############################################################
# 22. Write final Gene ID list
#
# EXACTLY:
#
# hum_geneID    type
#
# No header
############################################################

write.table(
  final_gene_list,
  OUT_FINAL,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)



############################################################
# 23. DONE
############################################################

cat("\n")
cat("========================================\n")
cat("ALL QC PASSED\n")
cat("========================================\n")


cat(
  "Input genes:",
  nrow(genes),
  "\n"
)


cat(
  "HGNC mapping:",
  sum(
    gene_map$HGNC_QC ==
      "PASS"
  ),
  "/",
  nrow(gene_map),
  "PASS\n"
)


cat(
  "NCBI validation:",
  sum(
    gene_map$NCBI_QC ==
      "PASS"
  ),
  "/",
  nrow(gene_map),
  "PASS\n"
)


cat(
  "Final Gene IDs:",
  nrow(final_gene_list),
  "\n"
)


cat("\nOutputs:\n")


cat(
  "  QC table:    ",
  OUT_QC,
  "\n"
)


cat(
  "  Problem list:",
  OUT_PROBLEM,
  "\n"
)


cat(
  "  Gene ID list:",
  OUT_FINAL,
  "\n"
)

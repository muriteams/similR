dat <- readLines("src/similarity.cpp")

line_ids <- which(grepl("^double (s|d)", dat[-1]) & grepl("^//' @alia", dat[-length(dat)]))

fun_keys  <- dat[grepl("^//'\\s+[-].+[(][0-9]+[)]", dat, perl = TRUE)]
fun_name  <- gsub("^double\\s+|[(]", "", dat[line_ids + 1])
fun_alias <- dat[line_ids] 

# Getting all the aliases
fun_keys <- regmatches(fun_keys, gregexpr("`(\"[a-zA-Z0-9&]+\")`", fun_keys))
fun_keys <- lapply(fun_keys, gsub, pattern = "`", replace = "", fixed = TRUE)


cbind(alias = fun_alias, sapply(fun_keys, paste, collapse = ","), fun_name)

# Writing function
program <- sapply(lapply(fun_keys, sprintf, fmt = "(s == %s)"), paste, collapse = " | ")

program <- sprintf("   else if (%s) fun = &%s;", program, fun_name)
program[1] <- gsub("else ", "", program[1])
program <- c(program, "   else Rcpp::stop(\"The statistic '%s' is not defined.\", s);")

program <- sprintf(
  "void getmetric(std::string s, funcPtr & fun) {\n%s\n   return;\n}",
  paste(program, collapse = "\n")
  )

program <- paste(
  "#ifndef SIMILR_GETMETRIC",
  "#define SIMILR_GETMETRIC 1",
  "typedef double (*funcPtr)(const std::vector< double > & table, bool normalized);",
  paste(sprintf("double %s(const std::vector< double > & table, bool normalize);", fun_name), collapse = "\n"),
  program,
  "#endif",
  sep = "\n"
)

cat(program, file = "src/getmetric.h", sep = "\n")

# Creating dataset
fun_alias <- gsub(".+@aliases ", "", fun_alias)
statistics <- data.frame(
  description = fun_alias,
  type        = ifelse(grepl("^d", fun_name), "distance", "similarity"),
  alias       = sapply(fun_keys, "[[", 1),
  other_alias = sapply(fun_keys, function(f) {
    paste(gsub("\"", "", f[-1]), collapse = ", ")
  }),
  stringsAsFactors = FALSE
)

statistics$alias <- gsub("\"", "", statistics$alias)
statistics

usethis::use_data(statistics, overwrite = TRUE)

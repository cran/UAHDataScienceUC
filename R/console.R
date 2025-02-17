#' @importFrom cli console_width
#' @importFrom cli cli_text
console.log <- function(txt, ...) {
  width <- cli::console_width()

  # Split the text into lines
  tmp <- strsplit(txt, '[\r\n]')[[1]]

  # Process every line
  lines <- NULL
  for (line in tmp) {
    # Get the whitespaces at the beginning of the line
    white <- get_whitespace(line)

    # Remove trailing whitespace from the line
    line <- substring(line, nchar(white) + 1)
    if (is.na(line)) line <- ''

    # Expand the tabs in the whitespace
    white <- gsub('\t', '  ', white)

    # Split the line into several lines which will fit in the console
    parts <- strsplit(
      line,
      paste0("(?<=.{", max(width - nchar(white), 1), "})"),
      perl = TRUE
    )[[1]]

    # Add the whitespace to the beginning of each part
    parts <- paste0(white, parts)

    # Add all parts to the processed lines list
    lines <- c(lines, parts)
  }

  # Make sure to print something...
  if (length(lines) < 1) lines <- ""

  # Print all processed lines
  for (line in lines) {
    message(line, ...)
  }
}

hline <- function() {
  console.log(strrep('_', cli::console_width()))
}

get_whitespace <- function(txt) {
  # Find whitespace at the beginning of the string
  fst_match <- gregexpr("^[ \t]*", txt)[[1]]

  # Extract whitespace
  white.length <- attr(fst_match, "match.length")
  substring(txt, 1, white.length)
}

format_subscript <- function(x) {
  # Map numbers to subscript Unicode characters
  subscripts <- c('\u2080', '\u2081', '\u2082', '\u2083', '\u2084',
                  '\u2085', '\u2086', '\u2087', '\u2088', '\u2089')
  chars <- strsplit(as.character(x), '')[[1]]
  paste0(sapply(chars, function(d) subscripts[as.numeric(d) + 1]), collapse='')
}

# Function to display distance formulas
display_distance_formula <- function(distance_type) {
  switch(distance_type,
         "euclidean" = {
           cli::cli_text("Euclidean Distance Formula:")
           cli::cli_text("d(x,y) = \u221A(\u2211{i=1}\u207F (x\u1D62 - y\u1D62)\u00B2)")
         },
         "manhattan" = {
           cli::cli_text("Manhattan Distance Formula:")
           cli::cli_text("d(x,y) = \u2211{i=1}\u207F |x\u1D62 - y\u1D62|")
         },
         "canberra" = {
           cli::cli_text("Canberra Distance Formula:")
           cli::cli_text("d(x,y) = \u2211{i=1}\u207F |x\u1D62 - y\u1D62|/(|x\u1D62| + |y\u1D62|)")
         },
         "chebyshev" = {
           cli::cli_text("Chebyshev Distance Formula:")
           cli::cli_text("d(x,y) = max{i=1}\u207F |x\u1D62 - y\u1D62|")
         },
         cli::cli_alert_warning("Unknown distance type")
  )
}

# Function to explain distance calculation with example
explain_distance_calculation <- function(x, y, distance_type) {
  # Get the actual distance using proxy package
  d <- proxy::dist(matrix(x, ncol=length(x)),
                   matrix(y, ncol=length(y)),
                   method=distance_type)
  # Display the formula
  display_distance_formula(distance_type)
  # Show calculation steps
  cli::cli_h2("Calculation Steps")
  switch(distance_type,
         "euclidean" = {
           diffs <- (x - y)^2
           cli::cli_text("1. Calculate squared differences:")
           for(i in seq_along(x)) {
             cli::cli_text("   (x{i} - y{i})\u00B2 = ({x[i]} - {y[i]})\u00B2 = {diffs[i]}",
                           i=format_subscript(i))
           }
           cli::cli_text("2. Sum the squared differences: {sum(diffs)}")
           cli::cli_text("3. Take the square root: \u221A{sum(diffs)} = {sqrt(sum(diffs))}")
         },
         "manhattan" = {
           diffs <- abs(x - y)
           cli::cli_text("1. Calculate absolute differences:")
           for(i in seq_along(x)) {
             cli::cli_text("   |x{i} - y{i}| = |{x[i]} - {y[i]}| = {diffs[i]}",
                           i=format_subscript(i))
           }
           cli::cli_text("2. Sum the absolute differences: {sum(diffs)}")
         }
  )
  cli::cli_text("\nFinal {distance_type} distance: {d}")
}

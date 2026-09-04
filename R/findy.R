#' What objects of this class are available
#'
#' Generic class finder
#'
#' Finds objects of the specified class in the EcoTest and EcoTestData R packages
#'
#' @param classy A class of object (character string, e.g. 'Fleet')
#' @param msg Print messages?
#' @param abc Logical, whether to alphabetize the results. By default, the function returns results found in
#' the global environment, then core openMSE packages, and any additional packages in argument `package`.
#' @author T. Carruthers
#' @seealso \link{Can} \link{Cant} \link{avail}
#' @examples
#' findy("ETData")
#' @export
findy <- function(classy, msg=TRUE, abc = FALSE) {

  get_funcs <- function(package, classy , msg) {
    pkgs <- search()
    search_package <- paste0("package:",package)
    funs <- NULL

    if (search_package %in% pkgs) {
      if (msg)
        cli::cli_alert('Searching for objects of class {.val {classy}} in package {.val {package}}')
      funs <- ls(search_package)[vapply(ls(search_package),
                                        getclass,
                                        logical(1),
                                        classy = classy)]
    } else {
      stop('Package ', package, ' not loaded. Use `library(', package, ')`', call. = FALSE)
    }
    funs
  }

  getclass <- function(x, classy) {
    return(any(class(get(x)) == classy)) # inherits(get(x), classy) - this gives a problem since we now inherit Stock etc in OM
  }

  global_funs <- ls(envir = .GlobalEnv)[vapply(ls(envir = .GlobalEnv), getclass, logical(1), classy = classy)]
  temp <- global_funs
  ETfuns <- get_funcs('EcoTest', classy, msg)
  temp <- c(temp, ETfuns)
  ETDfuns <- get_funcs('EcoTestData', classy, msg)
  temp <- c(temp, ETDfuns)

  if (length(temp) < 1) stop("No objects of class '", classy, "' found", call. = FALSE)

  out <- unique(temp)
  if (abc) out <- sort(out)

  return(out)

}


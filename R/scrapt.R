
#' Start SCRAPT shiny server
#'
#' @param threads  Number of threads to launch SCRAPT
#' @return
#' @export
#'
#' @examples

scrapt=function(threads=1){
    options(tidyverse.quiet = TRUE)
    shiny::shinyApp(ui = ui(), server = server(threads=threads))
}


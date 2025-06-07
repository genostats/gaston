#' Autosomes and X, Y, MT chromosomes
#' 
#' @description
#' Test if a chromosome id corresponds to an autosome or to X, Y, MT chromosomes
#' 
#' @param chr  A vector of chromosome ids
#' 
#' 
#' @details
#' These functions work by comparing the ids given in parameters with
#' the options \code{gaston.autosomes}, \code{gaston.chr.x}, \code{gaston.chr.y},
#' \code{gaston.chr.mt}.
#' 
#' For example, \code{is.autosome(chr)} is a short cut for
#' \code{chr \%in\% getOption("gaston.autosomes")}.
#' 
#' 
#' 
#' @return
#' A logical vector.
#' 
#' 
#' @keywords  Chromosome
#' @export is.autosome
is.autosome <- function(chr) chr %in% getOption("gaston.autosomes")

#' @rdname is.autosome
#' @export
is.chr.x    <- function(chr) chr %in% getOption("gaston.chr.x")

#' @rdname is.autosome
#' @export
is.chr.y    <- function(chr) chr %in% getOption("gaston.chr.y")

#' @rdname is.autosome
#' @export
is.chr.mt   <- function(chr) chr %in% getOption("gaston.chr.mt")



utils::globalVariables(
  c("MIC", "CI_Lower", "CI_Upper", "Delta_MIC", "Ratio_MIC",
    "Group1", "Group2", "Estimate", "obs_mic","rep",
    ".data","pretty_label","x_simplified", "vars")
)

#' autoplot method for mic_solve objects
#'
#' Quick ggplot2 visualisations of the main outputs from
#' [mic_solve()].  Three panels are supported:
#' * **"mic"**   - forest plot of group-wise MIC estimates with
#'   asymmetric CIs.
#' * **"delta"** - forest plot of deltaMIC pairwise differences.
#' * **"ratio"** - forest plot of MIC ratios (log scale).
#'
#' @param x        Variable for x axis plotting
#' @param object   An object returned by [mic_solve()].
#' @param type     One of `"mic"`, `"delta"`, or `"ratio"`, `"DoD_delta"`, or `"DoD_ratio"`.
#' @param color_by Optional column name used to color and dodge replicate
#'   points. Default: first column in `newdata`.
#' @param dot_size Size of the dots in the dotplot. Default: `0.5`.
#' @param ...      Additional arguments passed to [ggplot2::ggplot()].
#'
#' @return A **ggplot** object.
#' @author Carson Stacy
#' @examples
#' if (requireNamespace("ordinal", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   df <- data.frame(score  = ordered(sample(0:4, 120, TRUE)),
#'                    conc   = runif(120, 0, 4),
#'                    strain = factor(sample(c("A","B"), 120, TRUE)))
#'   fit <- ordinal::clm(score ~ strain * log1p(conc), data = df)
#'   res <- mic_solve(fit, expand.grid(strain = levels(df$strain)),
#'                    conc_name = "conc")
#'   ggplot2::autoplot(res, type = "mic")
#' }
#'
#' @name autoplot.mic_solve
#' @method autoplot mic_solve
#' @exportS3Method ggplot2::autoplot mic_solve
#' @aliases plot.mic_solve
#' @importFrom ggplot2 autoplot
#' @importFrom ggplot2 ggplot aes geom_pointrange geom_errorbar theme_minimal labs coord_flip theme element_text scale_y_log10
#' @importFrom stats ave relevel setNames formula terms
#'
autoplot.mic_solve <- function(object,
                               type = c("mic", "delta", "ratio",
                                        "DoD_delta", "DoD_ratio"),
                               x = NULL, color_by = NULL,# shape_by = NULL,
                               dot_size=0.5,
                               ...) {
  type <- match.arg(type)
  stopifnot(inherits(object, "mic_solve"))
  pd <- ggplot2::position_dodge(width = 0.5)

  if (type == "mic") {
    tbl <- object$mic_estimates
    group_vars <- names(object$newdata)

    if (is.null(x)) x <- names(object$newdata)[1]

    ## pick a color/shape column -------------------------------------------
    if (is.null(color_by)) color_by <- names(object$newdata)[min(2, ncol(object$newdata))]
    aes_extra <- list()
    if (color_by %in% names(tbl))
      aes_extra$color <- factor(tbl[[color_by]],
                                levels = levels(object$data[[color_by]]))
      aes_extra$fill <- #tbl[[color_by]]
        factor(tbl[[color_by]],
               levels = levels(object$data[[color_by]]))
    # if (!is.null(shape_by) && shape_by %in% names(tbl))
    #   aes_extra$shape <- tbl[[shape_by]]

    ## observed MIC per replicate -------------------------------------------
    d      <- object$data                   # raw data used in clm()
    conc_col <- if (identical(all.equal(object$transform_fun , identity), TRUE))
      object$conc_name else object$trans_name
    group_vars <- names(object$newdata)     # strain, treatment,

    ## ------------------------------------------------------------------ ##
    ## 1.  guarantee a replicate column

    # need a replicate id; if none supplied, treat each row as its own rep
    rep_col <- "rep"
    if (!rep_col %in% names(d)) {
      grp <- interaction(d[c(group_vars,conc_col)])
      d[[rep_col]] <- ave(d[[conc_col]], grp,
                          FUN = seq_along)
    }

    ## ------------------------------------------------------------------ ##
    ## 2.  observed MIC per replicate
    d$grp_key <- interaction(d[c(group_vars,rep_col)], drop = TRUE)

    split_data <- split(d, d$grp_key)

    score_var = as.character(attr(terms(formula(d)), "variables")[[2]])

    selected_rows <- lapply(split_data, function(group_df) {
      # Filter for minimum score
      filtered_df <- group_df[group_df[,score_var] > min(group_df[,score_var]), ]
      # Arrange by conc_var
      arranged_df <- filtered_df[order(filtered_df[[conc_col]],decreasing = T), ]
      # Select the first row
      arranged_df[1, ]
    })

    obs_tbl <- do.call(rbind, selected_rows)

    obs_tbl[[object$conc_name]] <- object$inv_transform_fun(obs_tbl[[conc_col]])

    ## assign color of points based on whether they are observed
    obs_tbl[["censored"]] <- ifelse(obs_tbl[[object$conc_name]] <
                                      max(obs_tbl[[object$conc_name]]),
                                  "Observed", "Not Observed")

    obs_tbl_censored <- obs_tbl
    obs_tbl_censored[[object$conc_name]] <- ifelse(obs_tbl_censored[["censored"]] == "Not Observed",
                                        obs_tbl_censored[[object$conc_name]],
                                        -1)

    obs_tbl_observed <- obs_tbl
    obs_tbl_observed[[object$conc_name]] <- ifelse(obs_tbl_observed[["censored"]] == "Observed",
                                        obs_tbl_observed[[object$conc_name]],
                                        -1)

    # hide extra points:
    obs_tbl_observed$visible  <- ifelse(obs_tbl_observed[[object$conc_name]] < 0, 0, 1)
    obs_tbl_censored$visible  <- ifelse(obs_tbl_censored[[object$conc_name]] < 0, 0, 0.5)
    bin_w <- max(obs_tbl[[object$conc_name]]) / 30

    ## ------------------------------------------------------------------ ##
    ## 2.  Plot Results

    ggplot2::ggplot(tbl,
                    ggplot2::aes(x = .data[[x]],
                                 y = MIC, ymin = CI_Lower, ymax = CI_Upper,
                                 !!!aes_extra)) +
      ggplot2::geom_col(position = pd, width = 0.5, alpha=0.5, color = NA) +
      # ggplot2::geom_pointrange(position = pd, size = 0.25, alpha=1) +
      ggplot2::geom_errorbar(position = pd, width = 0.25) +
      ggplot2::geom_dotplot(
                            data = obs_tbl_observed,
                            ggplot2::aes(x = obs_tbl_observed[[x]],
                                         y = obs_tbl_observed[[object$conc_name]],
                                         # alpha = visible,
                                         color = obs_tbl_observed[[color_by]],
                                         fill = obs_tbl_observed[[color_by]],
                                         group = interaction(obs_tbl[[x]],
                                                             obs_tbl[[color_by]]),
                            ),inherit.aes = FALSE,
                            binaxis = "y",      # Stack dots along the y-axis
                            # alpha = obs_tbl_observed$visible,
                            color = "black",
                            stackdir = "center", # Stack direction
                            dotsize = dot_size,      # Size of the dots
                            binwidth = max(obs_tbl[object$conc_name]/30),
                            position = ggplot2::position_dodge(width = pd$width)
                          ) +
      # for just the censored points
      ggplot2::geom_dotplot(
        data = obs_tbl_censored,
        ggplot2::aes(x = obs_tbl_censored[[x]],
                     y = obs_tbl_censored[[object$conc_name]],
                     color = obs_tbl_censored[[color_by]],
                     fill = obs_tbl_censored[[color_by]],
                     group = interaction(obs_tbl[[x]],
                                         obs_tbl[[color_by]])
        ), inherit.aes = FALSE,
        binaxis = "y",      # Stack dots along the y-axis
        color = "black",
        binpositions = "all",
        alpha = 0.5,#obs_tbl_censored$visible,
        stackdir = "center", # Stack direction
        dotsize = dot_size,      # Size of the dots
        binwidth = max(obs_tbl[object$conc_name]/30),
        position = ggplot2::position_dodge(width = pd$width)
      ) +
      ggplot2::scale_alpha(range = c(0, 1), guide = "none") +
      # ggplot2::scale_color_manual(values = levels(object$model$model[[color_by]])) +
      # ggplot2::scale_fill_manual(values = levels(object$model$model[[color_by]])) +
      # coord_flip() +
      # ggplot2::coord_cartesian(xlim = c(0, NA)) +  # applied while MIC is still on y
      ggplot2::coord_flip(ylim = c(0,#max(obs_tbl[object$conc_name])
                                   NA)
                          ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Group", y = "MIC",
                    color = color_by, fill = color_by, #shape = shape_by,
                    title = "Modelled MIC* with observed MICs") +
      ggplot2::theme(legend.position = "bottom")

  } else if (type == "delta") {

    ## ----- prettier deltaMIC forest plot  ------------------------------------
    tbl <- object$delta_mic_results

    ## --------------------------------------------------------------- ##
    ## 1.  explode Group1 / Group2 into factor columns
    vars <- names(object$newdata)              # e.g. "strain", "treatment"

    if (is.null(color_by)) color_by <- names(object$newdata)[min(2, ncol(object$newdata))]

    explode_delta <- function(col) {
      mat <- do.call(rbind,
                     strsplit(as.character(tbl[[col]]), ":", fixed = TRUE))
      out <- as.data.frame(mat, stringsAsFactors = FALSE)
      names(out) <- vars
      out
    }

    left  <- explode_delta("Group1")
    right <- explode_delta("Group2")
    names(left)  <- paste0(vars, "_1")
    names(right) <- paste0(vars, "_2")
    tbl <- cbind(tbl, left, right)

    for ( i in seq_along(vars) ) {
      tbl[[vars[i]]] <- ifelse(tbl[[paste0(vars[i], "_1")]]==
                               tbl[[paste0(vars[i], "_2")]],
                               tbl[[paste0(vars[i], "_1")]],
                               paste0(tbl[[paste0(vars[i], "_2")]],
                                      " - ",
                                      tbl[[paste0(vars[i], "_1")]]
                                      )
      )
    }

    explode <- function(x) {
      m <- do.call(rbind, strsplit(as.character(x),  ":", fixed = TRUE))
      out <- as.data.frame(m, stringsAsFactors = FALSE)
      names(out) <- vars
      out
    }

    g1 <- explode(tbl$Group1)
    g2 <- explode(tbl$Group2)

    ## 2. which variables match between the two groups? -----------------------
    same <- g1 == g2                           # nrow x n_vars logical

    tbl$facet <- apply(same, 1, function(r) {
      if (all(r))           "all"                       # every variable matches
      else if (!any(r))     paste0("Effect of other")   # none match
      else                  paste0("Effect of ",paste(vars[!r], collapse = "+"))  # e.g. "strain" or "strain+treatment"
    })

    tbl$facet <- relevel(as.factor(tbl$facet), "Effect of other")  # put "all" first




    tbl$x_simplified <- vapply(
      seq_len(nrow(tbl)),
      FUN.VALUE = character(1),
      function(i) {
        dif <- vars[g1[i, ] != g2[i, ]]
        sam <- vars[g1[i, ] == g2[i, ]]
        if (length(sam) == 0) sam <- NA  # if g2 and g1 are empty
        if (length(dif) == 0) {                 # all covariates same
          "no-diff"
        } else if (length(dif) == 1) {          # exactly one differs
          if (is.na(sam)) {
            # if g2 and g1 are empty
            sprintf("%s: %s - %s",
                    dif,
                    g2[i, dif],
                    g1[i, dif])                   # baseline (Group1)
            } else {
            # if g2 and g1 are not empty
            sprintf("%s: %s\n%s: %s - %s",
                    sam,
                    g2[i, sam],                   # "new" level (Group2)
                    dif,
                    g2[i, dif],                   # Group2 level first
                    g1[i, dif])                   # baseline (Group1)

            }
        } else {                                # >1 covariate differ
          paste(
                paste(g1[i, vars], collapse = ":"),#"",
                paste(g2[i, vars], collapse = ":"),
          sep = " - ")
        }
      }
    )



    ## 3. choose aesthetics ---------------------------------------------------
    if (is.null(color_by)) color_by <- names(object$newdata)[min(2, ncol(object$newdata))]
    aes_extra <- list()
    if (color_by %in% names(g1)) {             # use left-hand group's value
      #tbl[[color_by]] <- g1[[color_by]]
      aes_extra$color <- tbl[[color_by]]
      aes_extra$fill   <- tbl[[color_by]]
    }

    pd <- ggplot2::position_dodge(width = 0.6)

    ggplot2::ggplot(
      tbl,
      ggplot2::aes(x = x_simplified,#interaction(Group2, Group1, sep = " - "),
                   y = Delta_MIC,
                   ymin = CI_Lower, ymax = CI_Upper,
                   !!!aes_extra)
    ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          color = "grey60") +
      ggplot2::geom_col(position = pd, width = 0.5,
                        alpha = 0.35, color = NA) +
      ggplot2::geom_errorbar(position = pd, width = 0.25) +
      ggplot2::geom_point(position = pd, size = 2) +
      ggplot2::coord_flip() +
      ggplot2::facet_grid(rows=vars(factor(tbl$facet, levels = rev(levels(tbl$facet)))),
                          scales = "free_y", space='free') +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = NULL, y = "delta-MIC (absolute)",
                    color = color_by, fill = color_by,# shape = shape_by,
                    title = "Pairwise MIC differences (delta-MIC)") +
      ggplot2::theme(
        strip.text.y = ggplot2::element_text(face = "bold"),
        legend.position = if (is.null(color_by)) #&& is.null(shape_by))
          "none" else "bottom"
      )


  } else if (type == "ratio") {                           # ratio plot
    tbl <- object$ratio_mic_results

    if (is.null(tbl)) {
      stop("No pairwise MIC ratios found. Check that the model has at least two levels for each factor.")
    }

    ## --------------------------------------------------------------- ##
    ## 1.  explode Group1 / Group2 into factor columns
    vars <- names(object$newdata)              # e.g. "strain", "treatment"
    if (is.null(color_by)) color_by <- names(object$newdata)[min(2, ncol(object$newdata))]
    explode_mic <- function(col) {
      mat <- do.call(rbind,
                     strsplit(as.character(tbl[[col]]), ":", fixed = TRUE))
      out <- as.data.frame(mat, stringsAsFactors = FALSE)
      names(out) <- vars
      out
    }

    left  <- explode_mic("Group1")
    right <- explode_mic("Group2")
    names(left)  <- paste0(vars, "_1")
    names(right) <- paste0(vars, "_2")
    tbl <- cbind(tbl, left, right)

    for ( i in seq_along(vars) ) {
      tbl[[vars[i]]] <- ifelse(tbl[[paste0(vars[i], "_1")]]==
                                 tbl[[paste0(vars[i], "_2")]],
                               tbl[[paste0(vars[i], "_1")]],
                               paste0(tbl[[paste0(vars[i], "_2")]],
                                      " - ",
                                      tbl[[paste0(vars[i], "_1")]]
                               )
      )
    }

    explode <- function(x) {
      m <- do.call(rbind, strsplit(as.character(x),  ":", fixed = TRUE))
      out <- as.data.frame(m, stringsAsFactors = FALSE)
      names(out) <- vars
      out
    }

    g1 <- explode(tbl$Group1)
    g2 <- explode(tbl$Group2)

    ## 2. which variables match between the two groups? -----------------------
    same <- g1 == g2                           # nrow x n_vars logical

    tbl$facet <- apply(same, 1, function(r) {
      if (all(r))           "all"                       # every variable matches
      else if (!any(r))     paste0("Effect of other")   # none match
      else                  paste0("Effect of ",paste(vars[!r], collapse = "+"))  # e.g. "strain" or "strain+treatment"
    })

    tbl$facet <- relevel(as.factor(tbl$facet), "Effect of other")  # put "all" first

    # tbl$x_simplified <- vapply(
    #   seq_len(nrow(tbl)),
    #   FUN.VALUE = character(1),
    #   function(i) {
    #     dif <- vars[g1[i, ] != g2[i, ]]
    #     sam <- vars[g1[i, ] == g2[i, ]]
    #     if (length(sam) == 0) sam <- NA  # if g2 and g1 are empty
    #     if (length(dif) == 0) {                 # all covariates same
    #       "no-diff"
    #     } else if (length(dif) == 1) {          # exactly one differs
    #       if (is.na(sam)) {
    #         # if g2 and g1 are empty
    #         sprintf("%s: %s - %s",
    #                 dif,
    #                 g2[i, dif],
    #                 g1[i, dif])                   # baseline (Group1)
    #       } else {
    #         # if g2 and g1 are not empty
    #         sprintf("%s: %s\n%s: %s - %s",
    #                 sam,
    #                 g2[i, sam],                   # "new" level (Group2)
    #                 dif,
    #                 g2[i, dif],                   # Group2 level first
    #                 g1[i, dif])                   # baseline (Group1)
    #       }
    #     } else {                                # >1 covariate differ
    #       paste(
    #         paste(g1[i, vars], collapse = ":"),#"",
    #         paste(g2[i, vars], collapse = ":"),
    #         sep = " - ")
    #     }
    #   }
    # )

    # attempted fix for simplified x:

    tbl$x_simplified <- sapply(
      seq_len(nrow(tbl)),
      function(i) {
        # Use the actual vars that correspond to your model factors
        valid_vars <- vars[vars %in% names(g1) & vars %in% names(g2)]

        if (length(valid_vars) == 0) {
          return("no-vars")
        }

        # Find which variables differ and which are the same
        dif <- character(0)
        sam <- character(0)

        for (var in valid_vars) {
          g1_val <- g1[i, var]
          g2_val <- g2[i, var]

          # Skip if either value is NA or empty
          if (is.na(g1_val) || is.na(g2_val) ||
              length(g1_val) == 0 || length(g2_val) == 0) {
            next
          }

          if (g1_val == g2_val) {
            sam <- c(sam, var)
          } else {
            dif <- c(dif, var)
          }
        }

        if (length(dif) == 0) {
          # All variables are the same
          "no-diff"
        } else if (length(dif) == 1) {
          # Exactly one variable differs
          g1_val <- g1[i, dif]
          g2_val <- g2[i, dif]

          if (length(sam) == 0) {
            # No variables are the same
            sprintf("%s: %s vs. %s", dif, g2_val, g1_val)
          } else {
            # Some variables are the same
            sam_vals <- sapply(sam, function(s) g2[i, s])
            sam_collapsed <- paste(sprintf("%s: %s", sam, sam_vals), collapse = ", ")
            sprintf("%s\n%s: %s vs. %s", sam_collapsed, dif, g2_val, g1_val)
          }
        } else {
          # Multiple variables differ
          g1_vals <- sapply(valid_vars, function(var) g1[i, var])
          g2_vals <- sapply(valid_vars, function(var) g2[i, var])

          paste(
            paste(g1_vals, collapse = ":"),
            paste(g2_vals, collapse = ":"),
            sep = " - "
          )
        }
      },
      USE.NAMES = FALSE
    )


    if (is.null(color_by))
      color_by <- names(object$newdata)[min(2, ncol(object$newdata))]
    ## 3. choose aesthetics ---------------------------------------------------
    aes_extra <- list()
    if (color_by %in% names(g1)) {             # use left-hand group's value
      #tbl[[color_by]] <- g1[[color_by]]
      aes_extra$color <- tbl[[color_by]]
      aes_extra$fill   <- tbl[[color_by]]
    }

    pd <- ggplot2::position_dodge(width = 0.6)

    ggplot2::ggplot(
      tbl,
      ggplot2::aes(x = x_simplified,#interaction(Group2, Group1, sep = " - "),
                   y = log2(Ratio_MIC),
                   ymin = log2(CI_Lower), ymax = log2(CI_Upper),
                   !!!aes_extra)
    ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          color = "grey60") +
      ggplot2::geom_col(position = pd, width = 0.5,
                        alpha = 0.35, color = NA) +
      ggplot2::geom_errorbar(position = pd, width = 0.25) +
      ggplot2::geom_point(position = pd, size = 2) +
      ggplot2::coord_flip() +
      ggplot2::facet_grid(rows=vars(factor(tbl$facet, levels = rev(levels(tbl$facet)))),
                          scales = "free_y", space='free') +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Contrast", y = "MIC ratio (log scale)",
                    fill = color_by, color = color_by,# shape = shape_by,
           title = "Pairwise log2 MIC*", ...)






  } else if (type == "DoD_ratio") {
    ## ------------------Function to get pairs--------------------------------
    tbl  <- object$dod_ratio_results          # usually one row
    vars <- names(object$newdata)             # e.g. c("strain","treatment")

    if (nrow(tbl) == 0L) {
      stop("No pairwise difference of difference ratios found. Check that the model has at least two factors.")
    }

    ## ---------------------------------------------------------------
    ## 3.  aesthetic options
    aes_extra <- list(
      fill = tbl$label,
      color = tbl$label
      #fill = "#3182BD", color = "#3182BD"
    )

    ## ---------------------------------------------------------------
    ## 4.  plot
    ggplot2::ggplot(
      tbl,
      ggplot2::aes(x = .data[["label"]],
                   y = log2(Estimate),
                   ymin = log2(CI_Lower),
                   ymax = log2(CI_Upper),
                   !!!aes_extra)
    ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          color = "grey60") +
      ggplot2::geom_errorbar(width = .25, linewidth = .6) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_col(width = 0.5, alpha = 0.35, color = NA) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = NULL,
                    y = "Difference-of-differences (log2 ratio)",
                    title = "DoD on ratio scale") +
      ggplot2::theme(
        axis.text.y  = ggplot2::element_text(face = "bold"),
        legend.position = "none",
        plot.margin  = ggplot2::margin(r = 20)
      )


  } else if (type == "DoD_delta") {

    tbl  <- object$dod_delta_results          # usually one row
    vars <- names(object$newdata)             # e.g. c("strain","treatment")

    if (is.null(tbl)) {
      stop("No pairwise difference of differences found. Check that the model has at least two factors.")
    }

    aes_extra <- list(
      fill = tbl$label,
      color = tbl$label
      #fill = "#3182BD", color = "#3182BD"
      )

    ####  plot
    ggplot2::ggplot(
      tbl,
      ggplot2::aes(x = .data[["label"]],
                   y = Estimate,
                   ymin = CI_Lower,
                   ymax = CI_Upper,
                   !!!aes_extra)
    ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          color = "grey60") +
      ggplot2::geom_errorbar(width = .25, linewidth = .6) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_col(width = 0.5, alpha = 0.35, color = NA) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = NULL,
                    y = "Difference-of-differences",
                    title = "DoD of raw MIC") +
      ggplot2::theme(
        axis.text.y  = ggplot2::element_text(face = "bold"),
        legend.position = "none",
        plot.margin  = ggplot2::margin(r = 20)
      )
  }
}


utils::globalVariables(
  c("MIC", "CI_Lower", "CI_Upper", "Delta_MIC", "Ratio_MIC",
    "Group1", "Group2", "Estimate", "obs_mic","rep",
    ".data","pretty_label","x_simplified")
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
#' @importFrom stats ave relevel setNames
#'
autoplot.mic_solve <- function(object,
                               type = c("mic", "delta", "ratio",
                                        "DoD_delta", "DoD_ratio"),
                               x = NULL, color_by = NULL,# shape_by = NULL,
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
      aes_extra$color <- tbl[[color_by]]
      aes_extra$fill <- tbl[[color_by]]
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

    selected_rows <- lapply(split_data, function(group_df) {
      # Filter for minimum score
      filtered_df <- group_df[group_df$score > min(group_df$score), ]
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
    ## ------------------------------------------------------------------ ##
    ## 2.  Plot Results

    ggplot2::ggplot(tbl,
                    ggplot2::aes(x = .data[[x]],
                                 y = MIC, ymin = CI_Lower, ymax = CI_Upper,
                                 !!!aes_extra)) +
      ggplot2::geom_col(position = pd, width = 0.5, alpha=0.5, color = NA) +
      ggplot2::geom_pointrange(position = pd, size = 0.25, alpha=1) +
      ggplot2::geom_dotplot(
                            data = obs_tbl_observed,
                            ggplot2::aes(x = obs_tbl_observed[[x]],
                                         y = obs_tbl_observed[[object$conc_name]],
                                         # color = obs_tbl_observed[[color_by]],
                                         fill = obs_tbl_observed[[color_by]],
                                         group = interaction(obs_tbl[[x]],
                                                             obs_tbl[[color_by]]),
                            ),inherit.aes = FALSE,
                            binaxis = "y",      # Stack dots along the y-axis
                            color = "black",
                            stackdir = "center", # Stack direction
                            dotsize = 0.5,      # Size of the dots
                            binwidth = max(obs_tbl[object$conc_name]/30),
                            position = ggplot2::position_dodge(width = pd$width)
                          ) +
      # for just the censored points
      ggplot2::geom_dotplot(
        data = obs_tbl_censored,
        ggplot2::aes(x = obs_tbl_censored[[x]],
                     y = obs_tbl_censored[[object$conc_name]],
                     #color = obs_tbl_censored[[color_by]],
                     fill = obs_tbl_censored[[color_by]],
                     group = interaction(obs_tbl[[x]],
                                         obs_tbl[[color_by]] )
        ), inherit.aes = FALSE,
        binaxis = "y",      # Stack dots along the y-axis
        color = "black",
        binpositions = "all",
        alpha = 0.5,
        stackdir = "center", # Stack direction
        dotsize = 0.5,      # Size of the dots
        binwidth = max(obs_tbl[object$conc_name]/30),
        position = ggplot2::position_dodge(width = pd$width)
      ) +
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

    explode_ <- function(col) {
      mat <- do.call(rbind,
                     strsplit(as.character(tbl[[col]]), ":", fixed = TRUE))
      out <- as.data.frame(mat, stringsAsFactors = FALSE)
      names(out) <- vars
      out
    }

    left  <- explode_("Group1")
    right <- explode_("Group2")
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
        if (length(dif) == 0) {                 # all covariates same
          "no-diff"
        } else if (length(dif) == 1) {          # exactly one differs
          sprintf("%s: %s\n%s: %s - %s",
                  sam,
                  g2[i, sam],                   # "new" level (Group2)
                  dif,
                  g2[i, dif],                   # Group2 level first
                  g1[i, dif])                   # baseline (Group1)
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
    # if (!is.null(shape_by) && shape_by %in% names(g1)) {
    #   #tbl[[shape_by]]  <- g1[[shape_by]]
    #   aes_extra$shape  <- tbl[[shape_by]]
    # }

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
      ggplot2::facet_wrap(~ factor(facet, levels = rev(levels(facet))),
                          scales = "free_y", dir = "v") +
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

    ## --------------------------------------------------------------- ##
    ## 1.  explode Group1 / Group2 into factor columns
    vars <- names(object$newdata)              # e.g. "strain", "treatment"
    if (is.null(color_by)) color_by <- names(object$newdata)[min(2, ncol(object$newdata))]
    explode <- function(col) {
      mat <- do.call(rbind,
                     strsplit(as.character(tbl[[col]]), ":", fixed = TRUE))
      out <- as.data.frame(mat, stringsAsFactors = FALSE)
      names(out) <- vars
      out
    }

    left  <- explode("Group1")
    right <- explode("Group2")
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
        if (length(dif) == 0) {                 # all covariates same
          "no-diff"
        } else if (length(dif) == 1) {          # exactly one differs
          sprintf("%s: %s\n%s: %s - %s",
                  sam,
                  g2[i, sam],                   # "new" level (Group2)
                  dif,
                  g2[i, dif],                   # Group2 level first
                  g1[i, dif])                   # baseline (Group1)
        } else {                                # >1 covariate differ
          paste(
            paste(g1[i, vars], collapse = ":"),#"",
            paste(g2[i, vars], collapse = ":"),
            sep = " - ")
        }
      }
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
    # if (!is.null(shape_by) && shape_by %in% names(g1)) {
    #   #tbl[[shape_by]]  <- g1[[shape_by]]
    #   aes_extra$shape  <- tbl[[shape_by]]
    # }

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
      ggplot2::facet_wrap(~ factor(facet, levels = rev(levels(facet))),
                          scales = "free_y", dir = "v") +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Contrast", y = "MIC ratio (log scale)",
                    fill = color_by, color = color_by,# shape = shape_by,
           title = "Pairwise log2 MIC*", ...)






  } else if (type == "DoD_ratio") {
    # dat <- object$dod_ratio_results
    # ggplot2::ggplot(dat, ggplot2::aes(x = "Change in\ Response",
    #                                   y = log2(Estimate), ymin = log2(CI_Lower), ymax = log2(CI_Upper))) +
    #   ggplot2::geom_pointrange() +
    #   #ggplot2::scale_y_log10() +
    #   ggplot2::theme_minimal() + ggplot2::coord_flip() +
    #   ggplot2::labs(x = "DoD FC", y = "DoD ratio (log scale)",
    #                 title = "Difference of Differences", ...)
    ## ------------------Function to get pairs--------------------------------

    tbl  <- object$dod_ratio_results          # usually one row
    vars <- names(object$newdata)             # e.g. c("strain","treatment")

    ## ---------------------------------------------------------------
    ## 1.  recover the four corner groups in the 2x2 factorial design
    g1 <- strsplit(object$delta_mic_results$Group1[1], ":", fixed = TRUE)[[1]]
    g2 <- strsplit(object$delta_mic_results$Group2[2], ":", fixed = TRUE)[[1]]
    # g1 = baseline (WT / None), g2 = other corner (Mut / Salt)

    baseline <- setNames(g1, vars)            # WT , None
    contrast <- setNames(g2, vars)            # Mut, Salt

    ## ---------------------------------------------------------------
    ## 2.  build a descriptive label
    label <- sprintf(
      "%s - %s vs. %s - %s",
      #vars[2],             # treatment (column 2 of newdata)
      contrast[vars[2]], baseline[vars[2]],
      # vars[1],             # strain (column 1 of newdata)
      contrast[vars[1]], baseline[vars[1]]
    )
    tbl$pretty_label <- label                  # one row, one label

    ## ---------------------------------------------------------------
    ## 3.  (existing color logic; keep as-is or simplify)
    aes_extra <- list(fill = "#3182BD", color = "#3182BD")

    ## ---------------------------------------------------------------
    ## 4.  plot
    ggplot2::ggplot(
      tbl,
      ggplot2::aes(x = pretty_label,
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
    # dat <- object$dod_delta_results
    # ggplot2::ggplot(dat, ggplot2::aes(x = "Change in\ Response",
    #                                   y = Estimate, ymin = CI_Lower, ymax = CI_Upper)) +
    #   ggplot2::geom_pointrange() +
    #   ggplot2::theme_minimal() + ggplot2::coord_flip() +
    #   ggplot2::labs(x = "DoD MIC*", y = "DoD delta (absolute)",
    #                 title = "Difference of Differences", ...)
    tbl  <- object$dod_delta_results          # usually one row
    vars <- names(object$newdata)             # e.g. c("strain","treatment")

    ## ---------------------------------------------------------------
    ## 1.  recover the four corner groups in the 2x2 factorial design
    g1 <- strsplit(object$delta_mic_results$Group1[1], ":", fixed = TRUE)[[1]]
    g2 <- strsplit(object$delta_mic_results$Group2[2], ":", fixed = TRUE)[[1]]
    # g1 = baseline (WT / None), g2 = other corner (Mut / Salt)

    baseline <- setNames(g1, vars)            # WT , None
    contrast <- setNames(g2, vars)            # Mut, Salt

    ## ---------------------------------------------------------------
    ## 2.  build a descriptive label
    label <- sprintf(
      "%s - %s vs. %s - %s",
      #vars[2],             # treatment (column 2 of newdata)
      contrast[vars[2]], baseline[vars[2]],
      # vars[1],             # strain (column 1 of newdata)
      contrast[vars[1]], baseline[vars[1]]
    )
    tbl$pretty_label <- label                  # one row, one label

    ## ---------------------------------------------------------------
    ## 3.  (existing color logic; keep as-is or simplify)
    aes_extra <- list(fill = "#3182BD", color = "#3182BD")

    ## ---------------------------------------------------------------
    ## 4.  plot
    ggplot2::ggplot(
      tbl,
      ggplot2::aes(x = pretty_label,
                   y = log2(exp(Estimate)),
                   ymin = log2(exp(CI_Lower)),
                   ymax = log2(exp(CI_Upper)),
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
  }
}

